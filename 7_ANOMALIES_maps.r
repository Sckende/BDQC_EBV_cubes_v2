# Obj: production carte d'anomalies en vue de comparer les sorties de SDMs


library(rstac)
library(sf)
library(terra)
library(dplyr)


# --> choix espece
# espece <- "bonasa_umbellus"
# espece <- "thuja_occidentalis"
# espece <- "melospiza_georgiana"
espece <- "setophaga_virens"

# --> choix methode
# methode <- "Maxent"
# methode <- "brt"
methode <- "randomForest"

# chargement des rasters
max1 <- rast(paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/TdB_bench_maps/maps/CROPPED_QC_", espece, "_", methode, "_Predictors_Bias_Spatial.tif"))
max2 <- rast(paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/TdB_bench_maps/maps/CROPPED_QC_", espece, "_", methode, "_Predictors_Bias_noSpatial.tif"))
max3 <- rast(paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/TdB_bench_maps/maps/CROPPED_QC_", espece, "_", methode, "_Predictors_noBias_Spatial.tif"))
max4 <- rast(paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/TdB_bench_maps/maps/CROPPED_QC_", espece, "_", methode, "_Predictors_noBias_noSpatial.tif"))
max5 <- rast(paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/TdB_bench_maps/maps/CROPPED_QC_", espece, "_", methode, "_noPredictors_Bias_Spatial.tif"))
max6 <- rast(paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/TdB_bench_maps/maps/CROPPED_QC_", espece, "_", methode, "_noPredictors_noBias_Spatial.tif"))

# ---------- #
# Methode 1 - soustraction des valeurs des rasters
# ---------- #
# --> methodes comparables rapidement = Maxent, RF, BRT car production de carte en proba de presence - soustraction des valeurs entre deux rasters - valeurs comprises entre -1 & 1
x11()
par(mfrow = c(6, 6))
ras_ls <- list(max1, max2, max3, max4, max5, max6)
df <- data.frame()

for (i in 1:6) {
    for (j in 1:6) {
        ras_dif <- ras_ls[[i]] - ras_ls[[j]]
        plot(ras_dif, range = c(-1, 1))
        vect <- c(i, j, mean(abs(values(ras_dif)), na.rm = T))

        df <- rbind(df, vect)
        print(paste0("Comparaison ", i, " et ", j, " done"))
    }
}

df
names(df) <- c("i", "j", "mean_diff")
barplot(df$mean_diff, names.arg = paste0(df$i, df$j))

# ---------- #
# Methode 2 - dismo::nicheOverlap()
# ---------- #
# --> Utilisation du package dismo et production d'un indice D et/ou I inclus entre 0 (no overlap) et 1 (total overlap)
library(dismo)
library(raster)

?nicheOverlap
nicheOverlap(raster(max1), raster(max2))

mat_I <- matrix(nrow = 6, ncol = 6)
mat_D <- matrix(nrow = 6, ncol = 6)

for (i in 1:6) {
    for (j in 1:6) {
        mat_I[i, j] <- nicheOverlap(raster(ras_ls[[i]]), raster(ras_ls[[j]]), stat = "I")
        mat_D[i, j] <- nicheOverlap(raster(ras_ls[[i]]), raster(ras_ls[[j]]), stat = "D")
        print(paste0("Comparaison ", i, " et ", j, " done"))
    }
}

x11()
par(mfrow = c(6, 1))
for (i in 1:6) {
    barplot(mat_I[i, ])
}
