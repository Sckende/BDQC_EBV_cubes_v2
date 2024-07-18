# Obtention de tendance temporelle de richesse specifiques par division spatiales du qc

library(terra)
library(sf)
library(mapview)
library(RCurl)
library(rmapshaper)

# CERQ avec 5 nibveaux
# Ecozones
ecoz <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecozones.gpkg")
ecoz <- ms_simplify(ecoz)

# Ecoregions
ecor <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecoregions.gpkg")

# Ecoprovinces
# Ecodistricts
species <- read.table("INLA_species_list_final.txt", h = F)
 <- data.frame()
summ_rs <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("year", "poly_name", "poly_ID", "spe_rich")
colnames(summ_rs) <- x

for (year in 1992:2017) {
    path_list <- c()
    spe_rich <- 0
    for(i in 1:length(species[,1])){
        path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/", species[i, 1], "_range_", year, ".tif")
        if(url.exists(path)){
            path_list <- c(path_list, path)
        }
    stak <- rast(path_list)
    for(j in 1:nrow(ecoz)){
        test <- terra::crop(stak, ecoz[j,])
        sum(colSums(test[, -1]) > 0)
    }
    

    }

}

mapview(ecoz[7,])
