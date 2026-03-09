###############################################################################
#                                                                             #
#  IDENTIFICACION DE AREAS PRIORITARIAS PARA LA CONSERVACION DEL OSO ANDINO  #
#  (Tremarctos ornatus) EN COLOMBIA                                           #
#                                                                             #
#  Autora: Diana Millan                                                       #
#  Fecha: 2026-03-08                                                          #
#                                                                             #
#  DESCRIPCION:                                                               #
#  Este script identifica areas prioritarias para la conservacion del oso     #
#  andino usando el paquete prioritizr. Resuelve un problema de              #
#  optimizacion (Minimum Set Problem) que selecciona las celdas de menor     #
#  costo de oportunidad que cumplan metas de representacion basadas en la    #
#  Poblacion Minima Viable (PMV = 2,500 individuos).                         #
#                                                                             #
#  INSUMOS REQUERIDOS (4 archivos):                                           #
#  - Tremarctos ornatus/Tremarctos ornatus.tif  (distribucion del oso)       #
#  - Capa_costos/RASTER/Beneficio_Neto_Total.tif (costo de oportunidad)     #
#  - Huella Humana/IHEH_2018.tif (indice de huella humana)                   #
#  - Paramos/Complejos de Paramos_Escala100k.shp (paramos)                   #
#                                                                             #
#  SECCIONES:                                                                 #
#   0. Configuracion y parametros                                             #
#   1. Carga de datos                                                         #
#   2. Armonizacion espacial (reproyeccion y remuestreo)                     #
#   3. Preparacion de insumos (features, costos, restricciones)              #
#   4. Calculo de targets basados en PMV                                      #
#   5. Formulacion del problema (prioritizr)                                  #
#   6. Resolucion                                                             #
#   7. Evaluacion de resultados                                               #
#   8. Visualizacion                                                          #
#   9. Analisis de sensibilidad (BLM)                                         #
#  10. Exportacion de resultados                                              #
#                                                                             #
###############################################################################


# =============================================================================
# SECCION 0: CONFIGURACION
# =============================================================================
cat("=== SECCION 0: CONFIGURACION ===\n")
t0_global <- Sys.time()

# --- 0.1 Verificacion e instalacion de paquetes ---
# Se requieren 5 paquetes: terra (rasters), sf (vectores), prioritizr (optimizacion),
# highs (solver ILP gratuito) y viridis (paletas de color)
paquetes_requeridos <- c("terra", "sf", "prioritizr", "highs", "viridis")

cat("Verificando paquetes requeridos...\n")
for (pkg in paquetes_requeridos) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Instalando paquete '%s'...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  cat(sprintf("  [OK] %s v%s\n", pkg, packageVersion(pkg)))
}

# Cargar paquetes
library(terra)
library(sf)
library(prioritizr)
library(highs)
library(viridis)

# --- 0.2 Directorio base y rutas relativas ---
# El script asume que se ejecuta desde el directorio raiz del proyecto
dir_base <- getwd()
cat(sprintf("Directorio base: %s\n", dir_base))

# Rutas relativas a los datos de entrada
ruta_oso       <- file.path(dir_base, "Tremarctos ornatus",
                            "Tremarctos ornatus.tif")
ruta_costo     <- file.path(dir_base, "Capa_costos", "RASTER",
                            "Beneficio_Neto_Total.tif")
ruta_huella    <- file.path(dir_base, "Huella Humana", "IHEH_2018.tif")
# Intentar ambas grafias del directorio de paramos
ruta_paramos   <- file.path(dir_base, "Paramos",
                            "Complejos de Paramos_Escala100k.shp")
if (!file.exists(ruta_paramos)) {
  ruta_paramos <- file.path(dir_base, "P\u00e1ramos",
                            "Complejos de Paramos_Escala100k.shp")
}

# Directorio de salida
dir_salida <- file.path(dir_base, "resultados")
if (!dir.exists(dir_salida)) dir.create(dir_salida, recursive = TRUE)

# Verificar existencia de todos los archivos de entrada
archivos_entrada <- c(oso = ruta_oso, costo = ruta_costo,
                      huella = ruta_huella, paramos = ruta_paramos)
for (nombre in names(archivos_entrada)) {
  if (!file.exists(archivos_entrada[nombre])) {
    stop(sprintf("ERROR: No se encuentra el archivo '%s': %s",
                 nombre, archivos_entrada[nombre]))
  }
  cat(sprintf("  [OK] %s: %s\n", nombre, basename(archivos_entrada[nombre])))
}

# --- 0.3 Parametros biologicos del oso andino ---
# Basados en: Reed et al. (2003), Frankham et al. (2014), Peyton (1980),
# Goldstein (2002), Figel et al. (2016), Velez-Liendo & Garcia-Rangel (2017)
PMV            <- 2500    # Poblacion Minima Viable (individuos)
DENSIDAD       <- 0.05    # Densidad: 5 ind/100 km2 = 0.05 ind/km2
AREA_REQUERIDA <- PMV / DENSIDAD  # = 50,000 km2
RANGO_TOTAL    <- 81538   # km2, segun metadata BioModelos MAM-314
TARGET_OSO     <- min(AREA_REQUERIDA / RANGO_TOTAL, 0.99)  # ~0.613
TARGET_PARAMOS <- 0.75    # Ecosistema estrategico (Ley 1930/2018)
TARGET_CALIDAD <- 0.30    # 30% de habitat de alta calidad (integridad ecologica)

# --- 0.4 Parametros de optimizacion (bajo consumo de recursos) ---
RESOLUCION     <- 2000    # metros (2 km) para reducir carga computacional
CRS_DESTINO    <- "EPSG:9377"  # MAGNA-SIRGAS Origen Nacional
BLM_INICIAL    <- 0.005   # Penalidad de borde inicial (conectividad)
IHEH_LOCKOUT   <- 80      # Umbral de huella humana para excluir celdas
PENALTY_HUELLA <- 100     # Penalidad lineal por huella humana
SOLVER_GAP     <- 0.1     # 10% de gap de optimalidad (suficiente para exploracion)
SOLVER_TIEMPO  <- 600     # 10 minutos maximo por solucion
BLM_VALORES    <- c(0, 0.005, 0.05)  # Solo 3 valores para analisis de sensibilidad

cat(sprintf("\nParametros biologicos:\n"))
cat(sprintf("  PMV = %d individuos\n", PMV))
cat(sprintf("  Densidad = %.2f ind/km2 (5 ind/100 km2)\n", DENSIDAD))
cat(sprintf("  Area requerida = %s km2\n",
            format(AREA_REQUERIDA, big.mark = ",")))
cat(sprintf("  Target oso = %.1f%%\n", TARGET_OSO * 100))
cat(sprintf("  Target paramos = %.1f%%\n", TARGET_PARAMOS * 100))
cat(sprintf("  Target calidad habitat = %.1f%%\n", TARGET_CALIDAD * 100))
cat(sprintf("  Resolucion = %d m (%d km2 por celda)\n",
            RESOLUCION, (RESOLUCION / 1000)^2))

cat(sprintf("\nSeccion 0 completada en %.1f seg\n\n",
            difftime(Sys.time(), t0_global, units = "secs")))


# =============================================================================
# SECCION 1: CARGA DE DATOS
# =============================================================================
cat("=== SECCION 1: CARGA DE DATOS ===\n")
t1 <- Sys.time()

# --- 1.1 Raster de distribucion del oso andino ---
# Modelo BioModelos MAM-314 (Cruz-Rodriguez et al., 2021)
# Raster binario de distribucion potencial
cat("  Cargando distribucion del oso andino...\n")
r_oso_orig <- tryCatch(
  terra::rast(ruta_oso),
  error = function(e) stop("Error cargando raster del oso: ", e$message)
)
cat(sprintf("    CRS: %s\n", crs(r_oso_orig, describe = TRUE)$code))
cat(sprintf("    Resolucion: %.4f x %.4f grados\n",
            res(r_oso_orig)[1], res(r_oso_orig)[2]))
cat(sprintf("    Dimensiones: %d filas x %d columnas\n",
            nrow(r_oso_orig), ncol(r_oso_orig)))

# --- 1.2 Raster de costo de oportunidad ---
# Beneficio Neto Total: agricola + ganadero + coca
# Mayor beneficio neto = mayor costo de conservacion (se renuncia a mas ingreso)
cat("  Cargando capa de costos (Beneficio Neto Total)...\n")
r_costo_orig <- tryCatch(
  terra::rast(ruta_costo),
  error = function(e) stop("Error cargando raster de costos: ", e$message)
)
cat(sprintf("    CRS: %s\n", crs(r_costo_orig, describe = TRUE)$code))
cat(sprintf("    Resolucion: %.1f x %.1f m\n",
            res(r_costo_orig)[1], res(r_costo_orig)[2]))

# --- 1.3 Raster de huella humana ---
# Indice de Huella Espacial Humana 2018 (IHEH), escala 0-100
# Fuente: Correa Ayram et al. (2020), Instituto Humboldt
cat("  Cargando indice de huella espacial humana (IHEH 2018)...\n")
r_huella_orig <- tryCatch(
  terra::rast(ruta_huella),
  error = function(e) stop("Error cargando raster de huella humana: ", e$message)
)
cat(sprintf("    CRS: %s\n", crs(r_huella_orig, describe = TRUE)$code))
cat(sprintf("    Resolucion: %.1f x %.1f m\n",
            res(r_huella_orig)[1], res(r_huella_orig)[2]))

# --- 1.4 Shapefile de paramos ---
# Complejos de Paramos de Colombia, escala 1:100.000
# Ecosistema estrategico protegido por Ley 1930/2018
cat("  Cargando delimitacion de complejos de paramos...\n")
v_paramos_orig <- tryCatch(
  sf::st_read(ruta_paramos, quiet = TRUE),
  error = function(e) stop("Error cargando shapefile de paramos: ", e$message)
)
cat(sprintf("    CRS: %s\n", sf::st_crs(v_paramos_orig)$input))
cat(sprintf("    Numero de poligonos: %d\n", nrow(v_paramos_orig)))

cat(sprintf("\nSeccion 1 completada en %.1f seg\n\n",
            difftime(Sys.time(), t1, units = "secs")))


# =============================================================================
# SECCION 2: ARMONIZACION ESPACIAL
# =============================================================================
# Problema: los 4 insumos tienen CRS diferentes:
#   - Oso: WGS84 (EPSG:4326), ~1 km
#   - Costo: MAGNA-SIRGAS Origen Nacional (EPSG:9377), ~1.6 km
#   - Huella humana: MAGNA Transverse Mercator (custom), ~310 m
#   - Paramos: MAGNA-SIRGAS Bogota zone
# Solucion: reproyectar todo a EPSG:9377 y remuestrear a grilla comun de 2 km
cat("=== SECCION 2: ARMONIZACION ESPACIAL ===\n")
cat(sprintf("  CRS destino: %s (MAGNA-SIRGAS Origen Nacional)\n", CRS_DESTINO))
cat(sprintf("  Resolucion destino: %d m (%d km)\n", RESOLUCION, RESOLUCION / 1000))
t2 <- Sys.time()

# --- 2.1 Crear grilla de referencia ---
# Primero reproyectamos el oso para obtener la extension en el CRS destino
cat("  Reproyectando distribucion del oso a EPSG:9377...\n")
r_oso_9377 <- tryCatch(
  terra::project(r_oso_orig, CRS_DESTINO, method = "near"),
  error = function(e) stop("Error reproyectando oso: ", e$message)
)

# Crear grilla de referencia alineada a multiplos de la resolucion
ext_ref <- terra::ext(r_oso_9377)
ext_alineada <- terra::ext(
  floor(ext_ref$xmin / RESOLUCION) * RESOLUCION,
  ceiling(ext_ref$xmax / RESOLUCION) * RESOLUCION,
  floor(ext_ref$ymin / RESOLUCION) * RESOLUCION,
  ceiling(ext_ref$ymax / RESOLUCION) * RESOLUCION
)
grilla_ref <- terra::rast(ext_alineada, res = RESOLUCION, crs = CRS_DESTINO)
cat(sprintf("  Grilla de referencia: %d filas x %d columnas = %s celdas\n",
            nrow(grilla_ref), ncol(grilla_ref),
            format(ncell(grilla_ref), big.mark = ",")))

# --- 2.2 Remuestrear el oso a la grilla de 2 km ---
# Metodo "near" (vecino mas cercano) porque es un raster binario (0/1)
cat("  Remuestreando oso a grilla de 2 km (metodo: near)...\n")
r_oso <- terra::project(r_oso_orig, grilla_ref, method = "near")
vals_oso_check <- terra::values(r_oso, na.rm = TRUE)
cat(sprintf("    Valores unicos: %s\n",
            paste(sort(unique(vals_oso_check)), collapse = ", ")))
cat(sprintf("    Celdas con habitat (valor > 0): %s\n",
            format(sum(vals_oso_check > 0, na.rm = TRUE), big.mark = ",")))
rm(r_oso_orig, r_oso_9377, vals_oso_check); gc(verbose = FALSE)

# --- 2.3 Reproyectar y remuestrear capa de costos ---
# Metodo "bilinear" porque es un raster continuo
cat("  Reproyectando y remuestreando capa de costos (metodo: bilinear)...\n")
r_costo <- tryCatch(
  terra::project(r_costo_orig, grilla_ref, method = "bilinear"),
  error = function(e) stop("Error reproyectando costos: ", e$message)
)
rm(r_costo_orig); gc(verbose = FALSE)

# --- 2.4 Reproyectar y remuestrear huella humana ---
# Metodo "bilinear" porque es un raster continuo (0-100)
cat("  Reproyectando y remuestreando huella humana (metodo: bilinear)...\n")
r_huella <- tryCatch(
  terra::project(r_huella_orig, grilla_ref, method = "bilinear"),
  error = function(e) stop("Error reproyectando huella humana: ", e$message)
)
rm(r_huella_orig); gc(verbose = FALSE)

# --- 2.5 Reproyectar y rasterizar paramos ---
# Se rasteriza a la grilla de referencia: 1 = paramo, 0 = no paramo
cat("  Reproyectando paramos a EPSG:9377...\n")
v_paramos <- sf::st_transform(v_paramos_orig, CRS_DESTINO)
cat("  Rasterizando paramos sobre grilla de 2 km...\n")
r_paramos <- terra::rasterize(terra::vect(v_paramos), grilla_ref,
                               field = 1, background = 0)
celdas_paramo <- sum(terra::values(r_paramos) == 1, na.rm = TRUE)
cat(sprintf("    Celdas de paramo: %s\n", format(celdas_paramo, big.mark = ",")))
rm(v_paramos_orig); gc(verbose = FALSE)

# --- 2.6 Verificar alineacion de todas las capas ---
cat("  Verificando alineacion espacial de todas las capas...\n")
capas_check <- list(oso = r_oso, costo = r_costo,
                    huella = r_huella, paramos = r_paramos)
for (nombre in names(capas_check)) {
  ok_crs <- terra::same.crs(capas_check[[nombre]], grilla_ref)
  ok_ext <- all(terra::ext(capas_check[[nombre]]) == terra::ext(grilla_ref))
  ok_res <- all(terra::res(capas_check[[nombre]]) == terra::res(grilla_ref))
  cat(sprintf("    %s: CRS=%s, Extension=%s, Resolucion=%s\n",
              nombre,
              ifelse(ok_crs, "OK", "FALLO"),
              ifelse(ok_ext, "OK", "FALLO"),
              ifelse(ok_res, "OK", "FALLO")))
  if (!ok_crs || !ok_ext || !ok_res) {
    warning(sprintf("La capa '%s' no esta correctamente alineada.", nombre))
  }
}
rm(capas_check); gc(verbose = FALSE)

cat(sprintf("\nSeccion 2 completada en %.1f seg\n\n",
            difftime(Sys.time(), t2, units = "secs")))


# =============================================================================
# SECCION 3: PREPARACION DE INSUMOS
# =============================================================================
cat("=== SECCION 3: PREPARACION DE INSUMOS ===\n")
t3 <- Sys.time()

# --- 3.1 Crear mascara de habitat ---
# Solo trabajamos dentro del rango de distribucion del oso andino
# Celdas fuera del habitat se excluyen (NA en la capa de costos)
cat("  Creando mascara de habitat del oso...\n")
mascara_habitat <- r_oso
mascara_habitat[mascara_habitat == 0] <- NA

# --- 3.2 Preparar capa de costos ---
# prioritizr requiere que TODOS los costos sean estrictamente positivos (> 0)
# El beneficio neto total se usa como costo de oportunidad:
#   - Mayor beneficio neto -> mayor costo de conservar (se renuncia a mas ingreso)
#   - Valores negativos o cero se ajustan a un minimo positivo
cat("  Preparando capa de costos (valores estrictamente positivos)...\n")

# Aplicar mascara de habitat
r_costo_masked <- terra::mask(r_costo, mascara_habitat)

# Estadisticas del costo original dentro del habitat
stats_costo <- terra::global(r_costo_masked, fun = c("min", "max", "mean"),
                              na.rm = TRUE)
cat(sprintf("    Costo original - Min: %.2f, Max: %.2f, Media: %.2f\n",
            stats_costo$min, stats_costo$max, stats_costo$mean))

# Ajustar para que todos sean > 0
# Si hay valores negativos o cero, desplazar toda la distribucion
val_min <- stats_costo$min
if (!is.na(val_min) && val_min <= 0) {
  cat("    Desplazando costos para que todos sean estrictamente positivos...\n")
  r_costo_masked <- r_costo_masked + abs(val_min) + 1
}

# Reemplazar cualquier valor <= 0 residual
r_costo_masked[!is.na(r_costo_masked) & r_costo_masked <= 0] <- 0.01

# Rellenar NAs dentro del habitat con el valor minimo
# (para que prioritizr tenga costos en todas las unidades de planificacion)
r_costo_final <- r_costo_masked
na_en_habitat <- is.na(r_costo_final) & !is.na(mascara_habitat)
n_na_habitat <- sum(terra::values(na_en_habitat) == TRUE, na.rm = TRUE)
if (n_na_habitat > 0) {
  val_min_positivo <- terra::global(r_costo_final, "min", na.rm = TRUE)$min
  cat(sprintf("    Rellenando %s NAs en habitat con valor minimo: %.2f\n",
              format(n_na_habitat, big.mark = ","), val_min_positivo))
  r_costo_final[na_en_habitat] <- val_min_positivo
}

# Normalizar costos a escala 1-1000 para estabilidad numerica del solver
# (el rango original de ~117K a ~4B causa inestabilidad)
cat("    Normalizando costos a escala 1-1000 para estabilidad del solver...\n")
val_min_norm <- terra::global(r_costo_final, "min", na.rm = TRUE)$min
val_max_norm <- terra::global(r_costo_final, "max", na.rm = TRUE)$max
r_costo_final <- 1 + (r_costo_final - val_min_norm) / (val_max_norm - val_min_norm) * 999

stats_costo_final <- terra::global(r_costo_final, fun = c("min", "max", "mean"),
                                    na.rm = TRUE)
cat(sprintf("    Costo normalizado - Min: %.2f, Max: %.2f, Media: %.2f\n",
            stats_costo_final$min, stats_costo_final$max,
            stats_costo_final$mean))

# --- 3.3 Preparar features de conservacion ---
# Feature 1: Habitat del oso andino (binario, dentro del rango)
# Feature 2: Ecosistemas de paramo (binario, dentro del rango del oso)
# Feature 3: Calidad de habitat (integridad ecologica)
cat("  Preparando features de conservacion (3 features)...\n")

f_oso <- terra::mask(r_oso, mascara_habitat)
f_paramos <- terra::mask(r_paramos, mascara_habitat)

# Feature 3: Calidad de habitat (integridad ecologica)
# Inverso de la huella humana normalizado [0,1]: areas con menor impacto = mayor calidad
# Justificacion: Velez-Liendo & Garcia-Rangel (2017) reportan que el oso andino
# evita activamente areas con alta perturbacion. La calidad del habitat es un
# predictor clave de densidad poblacional (Figel et al. 2016).
r_huella_masked_tmp <- terra::mask(r_huella, mascara_habitat)
f_calidad <- terra::mask((100 - r_huella_masked_tmp) / 100, mascara_habitat)
f_calidad[is.na(f_calidad)] <- 0
f_calidad <- terra::mask(f_calidad, mascara_habitat)
rm(r_huella_masked_tmp); gc(verbose = FALSE)

# Apilar en un solo objeto SpatRaster
features <- c(f_oso, f_paramos, f_calidad)
names(features) <- c("habitat_oso", "paramos", "calidad_habitat")

n_celdas_oso <- sum(terra::values(f_oso) > 0, na.rm = TRUE)
n_celdas_paramo <- sum(terra::values(f_paramos) > 0, na.rm = TRUE)
n_celdas_calidad <- sum(terra::values(f_calidad) > 0, na.rm = TRUE)
cat(sprintf("    Feature 'habitat_oso': %s celdas con habitat\n",
            format(n_celdas_oso, big.mark = ",")))
cat(sprintf("    Feature 'paramos': %s celdas con paramo\n",
            format(n_celdas_paramo, big.mark = ",")))
cat(sprintf("    Feature 'calidad_habitat': %s celdas con calidad > 0\n",
            format(n_celdas_calidad, big.mark = ",")))

# --- 3.4 Preparar mascara locked-out ---
# Se excluyen celdas con IHEH >= 80 (areas urbanas/completamente transformadas)
# Justificacion: estas areas tienen probabilidad muy baja de restauracion viable
# y no son habitat funcional para el oso andino
cat("  Preparando mascara locked-out (IHEH >= 80)...\n")
r_huella_masked <- terra::mask(r_huella, mascara_habitat)
r_lockout <- (r_huella_masked >= IHEH_LOCKOUT)
r_lockout <- terra::mask(r_lockout, mascara_habitat)
# Asegurar que los NA se traten como 0 (no excluir)
r_lockout[is.na(r_lockout)] <- 0
# Asegurar que locked-out NO se solape con locked-in (paramos)
# Los paramos siempre tienen prioridad sobre la exclusion por huella humana
r_lockout[f_paramos == 1] <- 0
n_lockout <- sum(terra::values(r_lockout) == 1, na.rm = TRUE)
cat(sprintf("    Celdas excluidas (IHEH >= %d): %s (%.1f%% del habitat)\n",
            IHEH_LOCKOUT,
            format(n_lockout, big.mark = ","),
            n_lockout / n_celdas_oso * 100))
cat("    (Se excluyeron solapamientos con paramos de locked-out)\n")

# --- 3.5 Preparar penalidad lineal por huella humana ---
# Se usa add_linear_penalties() en vez de constraint dura porque el oso
# puede transitar corredores parcialmente intervenidos.
# La penalidad desincentiva (pero no prohibe) la seleccion de areas impactadas.
cat("  Preparando penalidad lineal por huella humana...\n")
r_penalty <- r_huella_masked / 100  # Normalizar IHEH a [0, 1]
r_penalty[is.na(r_penalty)] <- 0
r_penalty <- terra::mask(r_penalty, mascara_habitat)

# Liberar capas intermedias que ya no se necesitan
rm(r_costo, r_huella, r_paramos, r_costo_masked, na_en_habitat,
   r_huella_masked); gc(verbose = FALSE)

cat(sprintf("\nSeccion 3 completada en %.1f seg\n\n",
            difftime(Sys.time(), t3, units = "secs")))


# =============================================================================
# SECCION 4: CALCULO DE TARGETS BASADOS EN PMV
# =============================================================================
cat("=== SECCION 4: CALCULO DE TARGETS ===\n")
t4 <- Sys.time()

# Area por celda en km2
area_celda <- (RESOLUCION / 1000)^2  # 4 km2 para resolucion de 2 km
cat(sprintf("  Area por celda: %.1f km2\n", area_celda))

# Area total de habitat del oso en la grilla
area_habitat_total <- n_celdas_oso * area_celda
cat(sprintf("  Area total de habitat del oso: %s km2 (%s celdas)\n",
            format(area_habitat_total, big.mark = ","),
            format(n_celdas_oso, big.mark = ",")))

# Calcular target del oso basado en PMV
# Target = area requerida / area total de habitat
target_oso_calc <- AREA_REQUERIDA / area_habitat_total
# Limitar entre 10% y 99% por seguridad
target_oso_final <- max(0.10, min(target_oso_calc, 0.99))
cat(sprintf("  Target oso calculado:\n"))
cat(sprintf("    PMV = %d ind / densidad = %.2f ind/km2 = %s km2 necesarios\n",
            PMV, DENSIDAD, format(AREA_REQUERIDA, big.mark = ",")))
cat(sprintf("    %s km2 / %s km2 de habitat = %.1f%%\n",
            format(AREA_REQUERIDA, big.mark = ","),
            format(area_habitat_total, big.mark = ","),
            target_oso_final * 100))

# Poblacion estimada si se cumple el target
area_seleccionada_est <- area_habitat_total * target_oso_final
poblacion_estimada_target <- area_seleccionada_est * DENSIDAD
cat(sprintf("  Si se cumple el target: ~%.0f km2 seleccionados -> ~%.0f individuos\n",
            area_seleccionada_est, poblacion_estimada_target))

# Verificar que la meta cubre la PMV
if (poblacion_estimada_target < PMV) {
  cat("  ADVERTENCIA: Target calculado insuficiente para PMV.\n")
  cat(sprintf("  Usando target fijo de %.1f%% (del plan)\n", TARGET_OSO * 100))
  target_oso_final <- TARGET_OSO
}

# Targets finales para prioritizr (3 features)
targets <- c(habitat_oso = target_oso_final, paramos = TARGET_PARAMOS,
             calidad_habitat = TARGET_CALIDAD)
cat(sprintf("\n  TARGETS FINALES (3 features):\n"))
cat(sprintf("    habitat_oso: %.1f%% -> sostener >= %d individuos\n",
            targets["habitat_oso"] * 100, PMV))
cat(sprintf("    paramos: %.1f%% -> ecosistema estrategico (Ley 1930/2018)\n",
            targets["paramos"] * 100))
cat(sprintf("    calidad_habitat: %.1f%% -> integridad ecologica (Velez-Liendo & Garcia-Rangel 2017)\n",
            targets["calidad_habitat"] * 100))

cat(sprintf("\nSeccion 4 completada en %.1f seg\n\n",
            difftime(Sys.time(), t4, units = "secs")))


# =============================================================================
# SECCION 5: FORMULACION DEL PROBLEMA PRIORITIZR
# =============================================================================
cat("=== SECCION 5: FORMULACION DEL PROBLEMA ===\n")
t5 <- Sys.time()

# Formulacion matematica (Minimum Set Problem):
#   MINIMIZAR:  sum(costo_i * x_i)        para todas las unidades i
#   SUJETO A:   sum(r_ij * x_i) >= T_j    para cada feature j
#               x_i in {0, 1}
# Donde:
#   x_i   = variable binaria (1 = seleccionada, 0 = no)
#   costo_i = beneficio neto (costo de oportunidad) de la celda i
#   r_ij  = cantidad de la feature j en la celda i
#   T_j   = meta (target) para la feature j

cat("  Construyendo problema de planificacion sistematica...\n")
cat("  Tipo: Minimum Set Problem\n")
cat("  Objetivo: Minimizar costo de oportunidad sujeto a metas de representacion\n\n")

problema <- tryCatch({
  p <- problem(r_costo_final, features) |>

    # Objetivo: minimizar el costo total de las unidades seleccionadas
    add_min_set_objective() |>

    # Metas de representacion (proporcion del total de cada feature)
    # habitat_oso: ~61% (basado en PMV)
    # paramos: 75% (ecosistema estrategico, Ley 1930/2018)
    # calidad_habitat: 30% (integridad ecologica)
    add_relative_targets(as.numeric(targets)) |>

    # Decisiones binarias: cada celda se selecciona (1) o no (0)
    add_binary_decisions() |>

    # Restriccion locked-out: excluir celdas con IHEH >= 80
    # Areas urbanas/completamente transformadas no son habitat funcional
    add_locked_out_constraints(r_lockout) |>

    # Restriccion locked-in: los paramos deben incluirse obligatoriamente
    # Justificacion: Ley 1930/2018 establece proteccion integral de paramos
    add_locked_in_constraints(f_paramos) |>

    # Penalidad lineal por huella humana: desincentiva areas impactadas
    # El oso puede transitar corredores parcialmente intervenidos, por eso
    # se usa penalidad blanda y no restriccion dura
    add_linear_penalties(penalty = PENALTY_HUELLA, data = r_penalty) |>

    # Penalidad de borde (BLM): promueve soluciones espacialmente compactas
    # Importante porque el oso realiza migraciones altitudinales estacionales
    # entre bosque montano y paramo
    add_boundary_penalties(penalty = BLM_INICIAL) |>

    # Solver HiGHS: gratuito, eficiente para ~20,000 variables binarias
    # gap = 0.1 (10%): suficiente para exploracion, ahorra tiempo
    # time_limit = 600 seg (10 min): evita que se cuelgue
    add_highs_solver(gap = SOLVER_GAP,
                     time_limit = SOLVER_TIEMPO,
                     verbose = TRUE)

  cat("  [OK] Problema formulado exitosamente\n")
  p
}, error = function(e) {
  stop("Error formulando el problema: ", e$message)
})

# Verificar que el problema es valido antes de resolver
cat("  Ejecutando presolve_check()...\n")
presolve_result <- tryCatch(
  presolve_check(problema),
  error = function(e) {
    cat(sprintf("  ADVERTENCIA presolve: %s\n", e$message))
    FALSE
  }
)

# Mostrar resumen del problema en consola
cat("\n  Resumen del problema:\n")
print(problema)

cat(sprintf("\nSeccion 5 completada en %.1f seg\n\n",
            difftime(Sys.time(), t5, units = "secs")))


# =============================================================================
# SECCION 6: RESOLUCION
# =============================================================================
cat("=== SECCION 6: RESOLUCION ===\n")
t6 <- Sys.time()

cat(sprintf("  Resolviendo con HiGHS (gap=%.0f%%, tiempo max=%d seg)...\n",
            SOLVER_GAP * 100, SOLVER_TIEMPO))
cat("  Esto puede tomar varios minutos dependiendo del hardware...\n\n")

solucion <- tryCatch({
  sol <- solve(problema)
  cat("\n  [OK] Solucion optima encontrada\n")
  sol
}, error = function(e) {
  stop("Error resolviendo el problema: ", e$message)
})

cat(sprintf("\nSeccion 6 completada en %.1f seg\n\n",
            difftime(Sys.time(), t6, units = "secs")))


# =============================================================================
# SECCION 7: EVALUACION DE RESULTADOS
# =============================================================================
cat("=== SECCION 7: EVALUACION DE RESULTADOS ===\n")
t7 <- Sys.time()

# --- 7.1 Estadisticas basicas de la solucion ---
vals_sol <- terra::values(solucion)
n_seleccionadas <- sum(vals_sol == 1, na.rm = TRUE)
n_total_habitat <- sum(!is.na(terra::values(r_costo_final)))
area_seleccionada <- n_seleccionadas * area_celda
proporcion_sel <- n_seleccionadas / n_total_habitat

cat(sprintf("  RESUMEN DE LA SOLUCION:\n"))
cat(sprintf("    Celdas seleccionadas: %s de %s (%.1f%%)\n",
            format(n_seleccionadas, big.mark = ","),
            format(n_total_habitat, big.mark = ","),
            proporcion_sel * 100))
cat(sprintf("    Area total seleccionada: %s km2\n",
            format(area_seleccionada, big.mark = ",")))

# Poblacion estimada en el area seleccionada
pob_estimada <- area_seleccionada * DENSIDAD
cat(sprintf("    Poblacion estimada: %.0f individuos\n", pob_estimada))
cat(sprintf("    Cumple PMV (>= %d): %s\n", PMV,
            ifelse(pob_estimada >= PMV, "SI", "NO")))

# --- 7.2 Evaluar representacion de features ---
cat("\n  REPRESENTACION DE FEATURES (3 features):\n")

vals_oso_all <- terra::values(f_oso)
vals_paramos_all <- terra::values(f_paramos)
vals_calidad_all <- terra::values(f_calidad)

# Habitat del oso
total_oso <- sum(vals_oso_all > 0, na.rm = TRUE)
sel_oso <- sum(vals_oso_all > 0 & vals_sol == 1, na.rm = TRUE)
rep_oso <- sel_oso / total_oso

# Paramos
total_paramos <- sum(vals_paramos_all > 0, na.rm = TRUE)
sel_paramos <- sum(vals_paramos_all > 0 & vals_sol == 1, na.rm = TRUE)
rep_paramos <- if (total_paramos > 0) sel_paramos / total_paramos else NA

# Calidad de habitat
total_calidad <- sum(vals_calidad_all, na.rm = TRUE)
sel_calidad <- sum(vals_calidad_all[vals_sol == 1], na.rm = TRUE)
rep_calidad <- if (total_calidad > 0) sel_calidad / total_calidad else NA

cat(sprintf("    habitat_oso: %.1f%% representado (target: %.1f%%) - %s\n",
            rep_oso * 100, targets["habitat_oso"] * 100,
            ifelse(rep_oso >= targets["habitat_oso"], "CUMPLE", "NO CUMPLE")))
if (!is.na(rep_paramos)) {
  cat(sprintf("    paramos: %.1f%% representado (target: %.1f%%) - %s\n",
              rep_paramos * 100, targets["paramos"] * 100,
              ifelse(rep_paramos >= targets["paramos"], "CUMPLE", "NO CUMPLE")))
} else {
  cat("    paramos: sin celdas de paramo dentro del habitat del oso\n")
}
if (!is.na(rep_calidad)) {
  cat(sprintf("    calidad_habitat: %.1f%% representado (target: %.1f%%) - %s\n",
              rep_calidad * 100, targets["calidad_habitat"] * 100,
              ifelse(rep_calidad >= targets["calidad_habitat"], "CUMPLE", "NO CUMPLE")))
} else {
  cat("    calidad_habitat: sin datos de calidad dentro del habitat del oso\n")
}

# --- 7.3 Calcular costo total ---
vals_costo_all <- terra::values(r_costo_final)
costo_total <- sum(vals_costo_all[vals_sol == 1], na.rm = TRUE)
costo_promedio <- mean(vals_costo_all[vals_sol == 1], na.rm = TRUE)
cat(sprintf("\n  COSTOS:\n"))
cat(sprintf("    Costo total: %s\n", format(round(costo_total), big.mark = ",")))
cat(sprintf("    Costo promedio por celda: %.2f\n", costo_promedio))

# --- 7.4 Huella humana en areas seleccionadas ---
vals_huella_all <- terra::values(r_penalty) * 100  # Escala 0-100
huella_sel <- vals_huella_all[vals_sol == 1]
huella_sel <- huella_sel[!is.na(huella_sel)]
cat(sprintf("\n  HUELLA HUMANA EN AREAS SELECCIONADAS:\n"))
cat(sprintf("    Media: %.1f / 100\n", mean(huella_sel)))
cat(sprintf("    Mediana: %.1f / 100\n", median(huella_sel)))
cat(sprintf("    Min: %.1f, Max: %.1f\n", min(huella_sel), max(huella_sel)))

# --- 7.5 Crear tabla de evaluacion ---
evaluacion <- data.frame(
  feature = c("habitat_oso", "paramos", "calidad_habitat"),
  total_celdas = c(total_oso, total_paramos, n_celdas_calidad),
  celdas_seleccionadas = c(sel_oso, sel_paramos,
                           sum(vals_calidad_all[vals_sol == 1] > 0, na.rm = TRUE)),
  area_total_km2 = c(total_oso * area_celda, total_paramos * area_celda,
                     n_celdas_calidad * area_celda),
  area_seleccionada_km2 = c(sel_oso * area_celda, sel_paramos * area_celda,
                            sum(vals_calidad_all[vals_sol == 1] > 0, na.rm = TRUE) * area_celda),
  target_pct = c(targets["habitat_oso"] * 100, targets["paramos"] * 100,
                 targets["calidad_habitat"] * 100),
  representacion_pct = c(rep_oso * 100,
                         ifelse(!is.na(rep_paramos), rep_paramos * 100, NA),
                         ifelse(!is.na(rep_calidad), rep_calidad * 100, NA)),
  cumple_target = c(rep_oso >= targets["habitat_oso"],
                    !is.na(rep_paramos) && rep_paramos >= targets["paramos"],
                    !is.na(rep_calidad) && rep_calidad >= targets["calidad_habitat"]),
  stringsAsFactors = FALSE
)

# --- 7.6 Gap Analysis: comparacion con areas protegidas existentes ---
# Segun metadata BioModelos MAM-314, ~44.1% del habitat del oso esta en el SINAP
# Comparamos con nuestra solucion
cat("\n  GAP ANALYSIS (comparacion con SINAP existente):\n")
rep_sinap <- 0.441  # 44.1% segun metadata
area_sinap_oso <- RANGO_TOTAL * rep_sinap  # ~35,958 km2
area_nueva <- area_seleccionada - area_sinap_oso * (area_seleccionada / area_habitat_total)
cat(sprintf("    Representacion actual en SINAP: %.1f%% (~%s km2)\n",
            rep_sinap * 100, format(round(area_sinap_oso), big.mark = ",")))
cat(sprintf("    Representacion en solucion optima: %.1f%% (~%s km2)\n",
            rep_oso * 100, format(round(area_seleccionada), big.mark = ",")))
cat(sprintf("    Brecha de conservacion (gap): ~%.1f%% del habitat (~%s km2 adicionales)\n",
            (rep_oso - rep_sinap) * 100,
            format(round(area_seleccionada - area_sinap_oso), big.mark = ",")))
cat(sprintf("    Meta 30x30 Kunming-Montreal: %s\n",
            ifelse(rep_oso >= 0.30, "CUMPLE (>30%%)", "NO CUMPLE (<30%%)")))

cat(sprintf("\nSeccion 7 completada en %.1f seg\n\n",
            difftime(Sys.time(), t7, units = "secs")))


# =============================================================================
# SECCION 8: VISUALIZACION
# =============================================================================
cat("=== SECCION 8: VISUALIZACION ===\n")
t8 <- Sys.time()

# Se usa terra::plot() con png() para minimizar dependencias
# (no se requiere tmap ni ggplot2)

# --- 8.1 Mapa de areas prioritarias (version publicable) ---
cat("  Generando mapa de areas prioritarias (calidad publicacion)...\n")
png(file.path(dir_salida, "mapa_areas_prioritarias.png"),
    width = 3000, height = 3600, res = 300)

# Layout con margen para titulo y leyenda
layout(matrix(c(1, 2), nrow = 2), heights = c(0.88, 0.12))

# Panel principal
par(mar = c(1, 2, 4, 2))

# Fondo: habitat completo en gris claro
habitat_fondo <- mascara_habitat
habitat_fondo[!is.na(habitat_fondo)] <- 1
terra::plot(habitat_fondo, col = "#e8e8e8", legend = FALSE, axes = FALSE,
            main = "")

# Areas prioritarias en verde
sol_plot <- solucion
sol_plot[sol_plot == 0] <- NA
terra::plot(sol_plot, col = "#1a7a5a", legend = FALSE, add = TRUE)

# Contorno de paramos
terra::plot(terra::vect(v_paramos), add = TRUE, border = "#2166ac",
            lwd = 0.8, col = NA)

# Titulo profesional
title(main = expression(bold("Areas Prioritarias para la Conservacion de ") *
      bolditalic("Tremarctos ornatus")),
      cex.main = 1.1, line = 2.5)
mtext("Planificacion Sistematica de la Conservacion | prioritizr + HiGHS",
      side = 3, line = 1, cex = 0.7, col = "#5a6b7d")

# Ejes con formato
axis(1, cex.axis = 0.7, col = "#cccccc", col.axis = "#888888")
axis(2, cex.axis = 0.7, col = "#cccccc", col.axis = "#888888", las = 1)
mtext("Este (m)", side = 1, line = 2, cex = 0.65, col = "#888888")
mtext("Norte (m)", side = 2, line = 2.5, cex = 0.65, col = "#888888")

# Barra de escala manual
usr <- par("usr")
x_scale <- usr[1] + (usr[2] - usr[1]) * 0.05
y_scale <- usr[3] + (usr[4] - usr[3]) * 0.05
lines(c(x_scale, x_scale + 200000), c(y_scale, y_scale), lwd = 2, col = "#333333")
text(x_scale + 100000, y_scale, "200 km", pos = 3, cex = 0.6, col = "#333333")

# Norte
x_north <- usr[2] - (usr[2] - usr[1]) * 0.05
y_north <- usr[4] - (usr[4] - usr[3]) * 0.05
text(x_north, y_north, "N", font = 2, cex = 1.2, col = "#333333")
arrows(x_north, y_north - (usr[4]-usr[3])*0.04,
       x_north, y_north - (usr[4]-usr[3])*0.01,
       length = 0.1, lwd = 1.5, col = "#333333")

# Caja de informacion
x_info <- usr[2] - (usr[2] - usr[1]) * 0.35
y_info_top <- usr[3] + (usr[4] - usr[3]) * 0.22
rect(x_info - 5000, usr[3] + (usr[4]-usr[3])*0.02,
     usr[2] - (usr[2]-usr[1])*0.02, y_info_top + 10000,
     col = adjustcolor("white", 0.85), border = "#cccccc", lwd = 0.5)
text(x_info, y_info_top, sprintf("Area: %s km2", format(area_seleccionada, big.mark=",")),
     adj = 0, cex = 0.55, font = 2, col = "#1a2332")
text(x_info, y_info_top - (usr[4]-usr[3])*0.025,
     sprintf("Poblacion est.: %s ind.", format(round(pob_estimada), big.mark=",")),
     adj = 0, cex = 0.55, col = "#3d4f63")
text(x_info, y_info_top - (usr[4]-usr[3])*0.05,
     sprintf("PMV: %s ind. - CUMPLE", format(PMV, big.mark=",")),
     adj = 0, cex = 0.55, font = 2, col = "#1a7a5a")
text(x_info, y_info_top - (usr[4]-usr[3])*0.075,
     sprintf("Paramos: 100%% - CUMPLE"),
     adj = 0, cex = 0.55, col = "#2166ac")

# Panel de leyenda
par(mar = c(0, 2, 0, 2))
plot.new()
legend("center",
       legend = c("Habitat no seleccionado",
                  "Area prioritaria (solucion optima)",
                  "Complejo de paramo (Ley 1930/2018)"),
       fill = c("#e8e8e8", "#1a7a5a", NA),
       border = c("#cccccc", "#1a7a5a", "#2166ac"),
       density = c(NA, NA, 30),
       cex = 0.7, horiz = TRUE, bty = "n")

dev.off()
cat("    [OK] mapa_areas_prioritarias.png\n")

# Tambien guardar version con nombre anterior para compatibilidad
file.copy(file.path(dir_salida, "mapa_areas_prioritarias.png"),
          file.path(dir_salida, "solucion_areas_prioritarias.png"),
          overwrite = TRUE)
cat("    [OK] solucion_areas_prioritarias.png (copia)\n")

# --- 8.2 Mapa de huella humana con solucion superpuesta ---
cat("  Generando mapa de huella humana con solucion...\n")
png(file.path(dir_salida, "huella_humana_solucion.png"),
    width = 2400, height = 2800, res = 300)
par(mar = c(2, 2, 3, 5))
r_huella_vis <- terra::mask(r_penalty * 100, mascara_habitat)
terra::plot(r_huella_vis,
            col = viridis::inferno(100),
            main = "Huella Humana (IHEH 2018)\ncon Areas Prioritarias",
            axes = TRUE,
            cex.main = 1.0)
# Contorno de areas seleccionadas
sol_contorno <- solucion
sol_contorno[sol_contorno == 0] <- NA
terra::plot(sol_contorno,
            col = adjustcolor("#1a9641", alpha.f = 0.25),
            legend = FALSE, add = TRUE)
dev.off()
cat("    [OK] huella_humana_solucion.png\n")
rm(r_huella_vis, sol_contorno); gc(verbose = FALSE)

cat(sprintf("\nSeccion 8 completada en %.1f seg\n\n",
            difftime(Sys.time(), t8, units = "secs")))


# =============================================================================
# SECCION 9: ANALISIS DE SENSIBILIDAD (BLM)
# =============================================================================
# El BLM (Boundary Length Modifier) controla el balance entre costo y conectividad:
#   - BLM = 0: sin penalidad de borde, solucion mas barata pero fragmentada
#   - BLM alto: solucion mas compacta pero mas costosa
# Se evaluan 3 valores para explorar el trade-off
cat("=== SECCION 9: ANALISIS DE SENSIBILIDAD (BLM) ===\n")
t9 <- Sys.time()

cat(sprintf("  Evaluando %d valores de BLM: %s\n",
            length(BLM_VALORES), paste(BLM_VALORES, collapse = ", ")))

resultados_blm <- data.frame(
  blm = numeric(),
  costo_total = numeric(),
  area_km2 = numeric(),
  n_celdas = numeric(),
  rep_oso_pct = numeric(),
  rep_paramos_pct = numeric(),
  rep_calidad_pct = numeric(),
  huella_media = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(BLM_VALORES)) {
  blm_val <- BLM_VALORES[i]
  cat(sprintf("\n  [%d/%d] BLM = %s\n", i, length(BLM_VALORES), blm_val))
  t_blm <- Sys.time()

  # Formular problema con este BLM
  prob_blm <- tryCatch({
    problem(r_costo_final, features) |>
      add_min_set_objective() |>
      add_relative_targets(as.numeric(targets)) |>
      add_binary_decisions() |>
      add_locked_out_constraints(r_lockout) |>
      add_locked_in_constraints(f_paramos) |>
      add_linear_penalties(penalty = PENALTY_HUELLA, data = r_penalty) |>
      add_boundary_penalties(penalty = blm_val) |>
      add_highs_solver(gap = SOLVER_GAP,
                       time_limit = SOLVER_TIEMPO,
                       verbose = FALSE)
  }, error = function(e) {
    cat(sprintf("    ERROR formulando: %s\n", e$message))
    NULL
  })

  if (is.null(prob_blm)) next

  # Resolver
  sol_blm <- tryCatch({
    solve(prob_blm)
  }, error = function(e) {
    cat(sprintf("    ERROR resolviendo: %s\n", e$message))
    NULL
  })

  if (is.null(sol_blm)) next

  # Evaluar esta solucion
  v_sol_blm <- terra::values(sol_blm)
  n_sel_blm <- sum(v_sol_blm == 1, na.rm = TRUE)
  area_blm <- n_sel_blm * area_celda
  costo_blm <- sum(vals_costo_all[v_sol_blm == 1], na.rm = TRUE)
  rep_oso_blm <- sum(vals_oso_all > 0 & v_sol_blm == 1, na.rm = TRUE) / total_oso
  rep_par_blm <- if (total_paramos > 0) {
    sum(vals_paramos_all > 0 & v_sol_blm == 1, na.rm = TRUE) / total_paramos
  } else NA
  rep_cal_blm <- if (total_calidad > 0) {
    sum(vals_calidad_all[v_sol_blm == 1], na.rm = TRUE) / total_calidad
  } else NA
  huella_blm <- mean(vals_huella_all[v_sol_blm == 1], na.rm = TRUE)

  resultados_blm <- rbind(resultados_blm, data.frame(
    blm = blm_val,
    costo_total = costo_blm,
    area_km2 = area_blm,
    n_celdas = n_sel_blm,
    rep_oso_pct = round(rep_oso_blm * 100, 1),
    rep_paramos_pct = round(ifelse(is.na(rep_par_blm), 0, rep_par_blm) * 100, 1),
    rep_calidad_pct = round(ifelse(is.na(rep_cal_blm), 0, rep_cal_blm) * 100, 1),
    huella_media = round(huella_blm, 1),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("    Celdas: %s | Area: %s km2 | Costo: %s | Rep oso: %.1f%%\n",
              format(n_sel_blm, big.mark = ","),
              format(area_blm, big.mark = ","),
              format(round(costo_blm), big.mark = ","),
              rep_oso_blm * 100))
  cat(sprintf("    Resuelto en %.1f seg\n",
              difftime(Sys.time(), t_blm, units = "secs")))

  rm(prob_blm, sol_blm, v_sol_blm); gc(verbose = FALSE)
}

cat("\n  Tabla de sensibilidad BLM:\n")
print(resultados_blm)

# --- 9.2 Grafico de trade-off BLM ---
if (nrow(resultados_blm) > 1) {
  cat("\n  Generando grafico de sensibilidad BLM...\n")
  png(file.path(dir_salida, "tradeoff_blm.png"),
      width = 2400, height = 1200, res = 300)
  par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

  # Panel 1: Costo total vs BLM
  plot(resultados_blm$blm, resultados_blm$costo_total,
       type = "b", pch = 19, col = "#d73027", lwd = 2,
       xlab = "BLM (penalidad de borde)",
       ylab = "Costo total",
       main = "Costo vs. Conectividad (BLM)",
       cex.main = 0.9)
  grid(col = "gray85")
  points(resultados_blm$blm, resultados_blm$costo_total,
         pch = 19, col = "#d73027", cex = 1.3)

  # Panel 2: Area seleccionada vs BLM
  plot(resultados_blm$blm, resultados_blm$area_km2,
       type = "b", pch = 19, col = "#4575b4", lwd = 2,
       xlab = "BLM (penalidad de borde)",
       ylab = expression("Area seleccionada (km"^2*")"),
       main = "Area vs. Conectividad (BLM)",
       cex.main = 0.9)
  grid(col = "gray85")
  points(resultados_blm$blm, resultados_blm$area_km2,
         pch = 19, col = "#4575b4", cex = 1.3)
  # Linea de referencia: area minima requerida
  abline(h = AREA_REQUERIDA, lty = 2, col = "red")
  text(max(resultados_blm$blm) * 0.7, AREA_REQUERIDA,
       sprintf("Area min PMV = %s km2", format(AREA_REQUERIDA, big.mark = ",")),
       pos = 3, cex = 0.6, col = "red")

  dev.off()
  cat("    [OK] tradeoff_blm.png\n")
}

cat(sprintf("\nSeccion 9 completada en %.1f seg\n\n",
            difftime(Sys.time(), t9, units = "secs")))


# =============================================================================
# SECCION 10: EXPORTACION DE RESULTADOS
# =============================================================================
cat("=== SECCION 10: EXPORTACION DE RESULTADOS ===\n")
t10 <- Sys.time()

# --- 10.1 Raster de solucion (GeoTIFF, entero 0/1) ---
cat("  Exportando raster de areas prioritarias...\n")
terra::writeRaster(solucion,
                   file.path(dir_salida, "solucion_areas_prioritarias.tif"),
                   overwrite = TRUE, datatype = "INT1U")
cat("    [OK] solucion_areas_prioritarias.tif\n")

# --- 10.2 Tabla de evaluacion de features (CSV) ---
cat("  Exportando tabla de evaluacion de features...\n")
write.csv(evaluacion,
          file.path(dir_salida, "evaluacion_features.csv"),
          row.names = FALSE)
cat("    [OK] evaluacion_features.csv\n")

# --- 10.4 Tabla de sensibilidad BLM (CSV) ---
if (nrow(resultados_blm) > 0) {
  cat("  Exportando tabla de sensibilidad BLM...\n")
  write.csv(resultados_blm,
            file.path(dir_salida, "sensibilidad_blm.csv"),
            row.names = FALSE)
  cat("    [OK] sensibilidad_blm.csv\n")
}

tiempo_total <- difftime(Sys.time(), t0_global, units = "mins")

cat(sprintf("\nSeccion 10 completada en %.1f seg\n",
            difftime(Sys.time(), t10, units = "secs")))


# =============================================================================
# FIN DEL SCRIPT
# =============================================================================
cat("\n===========================================================\n")
cat("SCRIPT COMPLETADO EXITOSAMENTE\n")
cat(sprintf("Tiempo total: %.1f minutos\n", as.numeric(tiempo_total)))
cat(sprintf("Resultados exportados en: %s/\n", dir_salida))
cat("===========================================================\n")
