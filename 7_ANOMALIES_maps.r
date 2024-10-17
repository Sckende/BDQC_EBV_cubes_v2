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
# *** Methode 2 *** - dismo::nicheOverlap()
# ---------- #
# --> Utilisation du package dismo et production d'un indice D et/ou I inclus entre 0 (no overlap) et 1 (total overlap)
library(dismo)
library(raster)
library(ggcorplot)

?nicheOverlap
nicheOverlap(raster(max1), raster(max2))

ras_ls <- list(max1, max2, max3, max4, max5, max6)

mat_I <- matrix(nrow = 6, ncol = 6)
mat_D <- matrix(nrow = 6, ncol = 6)

for (i in 1:6) {
    for (j in 1:6) {
        mat_I[i, j] <- nicheOverlap(raster(ras_ls[[i]]), raster(ras_ls[[j]]), stat = "I")
        mat_D[i, j] <- nicheOverlap(raster(ras_ls[[i]]), raster(ras_ls[[j]]), stat = "D")
    }
}

colnames(mat_D) <- c("Pred_Bias_Spat", "Pred_Bias_noSpat", "Pred_noBias_Spat", "Pred_noBias_noSpat", "noPred_Bias_Spat", "noPred_noBias_Spat")
rownames(mat_D) <- c("Pred_Bias_Spat", "Pred_Bias_noSpat", "Pred_noBias_Spat", "Pred_noBias_noSpat", "noPred_Bias_Spat", "noPred_noBias_Spat")

# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
}
reorder_cormat <- function(cormat) {
    # Use correlation between variables as distance
    dd <- as.dist((1 - cormat) / 2)
    hc <- hclust(dd)
    cormat <- cormat[hc$order, hc$order]
}
# Visualisation
library(reshape2)
library(ggplot2)
library(corrplot) # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

low_mat <- get_lower_tri(mat_D)
low_mat2 <- melt(low_mat, na.rm = T)

# matrice reordered
mat_reord <- reorder_cormat(mat_D)
mat_reord2 <- get_upper_tri(mat_reord)
mat_reord3 <- melt(mat_reord2, na.rm = T)

ggplot(data = mat_reord3, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
        low = "#ffffff", high = "darkgreen", limit = c(0, 1), space = "Lab",
        name = "D index"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
        angle = 45, vjust = 1,
        size = 12, hjust = 1
    )) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4)

corrplot(mat_D, type = "lower", col = COL1("YlGn", n = 200), col.lim = c(0, 1), addgrid.col = "white", addCoef.col = "white", method = "circle", is.corr = F, title = paste0(espece, " - ", methode), diag = F, outline = F, order = "AOE", tl.col = "black")

# ---------- #
# Methode 3 - Wilson, 2011 - https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2011.00115.x#b73
# ---------- #
# --> calcul de la similarite/dissimilarite sur la base des distances euclidienne et d'Hellinger puis PCoA

# library(philentropy)
library(labdsv)

mat1 <- as.matrix(max1)
mat2 <- as.matrix(max2)

# code from supplementary information - Wilson 2011
HellingerDist <- function(mat1, mat2) {
    p1 <- sum(mat1, na.rm = T)
    p2 <- sum(mat2, na.rm = T)
    return(sqrt(0.5 * sum((sqrt(mat1 / p1) - sqrt(mat2 / p2))^2, na.rm = T)))
}
HellingerDist(mat1, mat2)

EuclideanDist <- function(mat1, mat2) {
    return(sqrt(sum((mat1 - mat2)^2, na.rm = T)))
}
EuclideanDist(mat1, mat2)

NormEuclideanDist <- function(mat1, mat2) {
    p1 <- mat1 / sum(mat1, na.rm = T)
    p2 <- mat2 / sum(mat2, na.rm = T)
    return(sqrt(sum((p1 - p2)^2, na.rm = T)))
}
NormEuclideanDist(mat1, mat2)

hell_mat <- matrix(nrow = 6, ncol = 6)
eucli_mat <- matrix(nrow = 6, ncol = 6)

for (i in 1:6) {
    for (j in 1:6) {
        hell_mat[i, j] <- HellingerDist(as.matrix(ras_ls[[i]]), as.matrix(ras_ls[[j]]))
        eucli_mat[i, j] <- EuclideanDist(as.matrix(ras_ls[[i]]), as.matrix(ras_ls[[j]]))
    }
}

hell_pco <- pco(as.dist(hell_mat))
eucli_pco <- pco(as.dist(eucli_mat))


# visualisation
pcoPlotFunc <- function(
    xpts, ypts, numTests, numSteps, numReps, aTitle, PlotChars, PlotCol,
    PlotTag) {
    plot(xpts, ypts, xlab = "Axis 1", ylab = "Axis 2", pch = PlotChars, col = PlotCol, main = aTitle, cex = 1.5)
    mtext(PlotTag, line = 1, font = 2, adj = 0)
    text(xpts[1], ypts[1], "Start", pos = 2, offset = 1)
    text(xpts[length(xpts)], ypts[length(xpts)], "End 1,3", pos = 4, offset = 1)
    n <- 1 + (numSteps - 1) + (numSteps * numReps)
    text(xpts[n], ypts[n], "End 2", pos = 4, offset = 1)
    j <- 1:numSteps
    for (i in 1:numTests)
    {
        k <- numSteps * (i - 1) + j
        lines(xpts[k], ypts[k], lty = 2)
    }
}
# Number of moving patterns or “sequences”
numSeqs <- 1
# Number of steps along a sequence
numSteps <- 1
# Number of replicates of a sequence to be generated
numReps <- 1

plotChr <- c(rep(16, numSteps * numReps), rep(17, numSteps * numReps), rep(3, numSteps * numReps))
c1 <- rgb(102, 102, 51, maxColorValue = 255)
c2 <- rgb(102, 102, 255, maxColorValue = 255)
c3 <- rgb(204, 204, 102, maxColorValue = 255)
plotCol <- rep(c(c1, c2, c3), numSeqs * numReps)

x11()
par(mfrow = c(1, 2))
# PCoA of Hellinger Distance
pcoPlotFunc(
    hell_pco$points[, 1], hell_pco$points[, 2], numSeqs * numReps, numSteps, numReps, "",
    plotChr, plotCol, "a. Hellinger"
)
# PCoA of Euclidean Distance
pcoPlotFunc(
    eucli_pco$points[, 1], eucli_pco$points[, 2], numSeqs * numReps, numSteps, numReps, "",
    plotChr, plotCol, "b. Euclidean"
)

# ==> to continue

# ---------- #
# Methode 4 - comparaison des histogrammes de valeurs de pixels
# ---------- #
x11()
par(mfrow = c(3, 2))
for (i in 1:6) {
    hist(values(ras_ls[[i]]), breaks = seq(0, 1, 0.01))
}
hist(values(max1))

#### PERMANOVA ####
library(vegan)

mat <- as.dist(hell_mat)

# ----> see the adonis2() function
