#### Prod carte RS from eBird maps (birds only)
# ----------------------------------------- #
# la liste d'especes correspond aux oiseaux nicheurs du Qc
# identique a celle utilisee par V. Bellavance
# une seule carte de RS produite car pas de tendance annuelle pour eBird

library(terra)
library(gdalcubes)
library(RCurl) # for check if url exists - url.exists()

# rast1 <- rast("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/Boulanger_map_converted/ABIE.BAL_Boulanger_EPSG_32198.tif")
# rast2 <- rast("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/Boulanger_map_converted/aggregate_ABIE.BAL_Boulanger_EPSG_32198.tif")

species <- read.csv("Boulanger_species_list_final.txt", header = F, col.names = c("species"))

rast_tot <- c()
for (i in 14:length(species$species)) {
    path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/Boulanger_map_converted/aggregate_", species[i, 1], "_Boulanger_EPSG_32198.tif")

    if (url.exists(path)) {
        ras <- rast(path)
        ras_agg <- aggregate(ras, fact = 2, fun = "sum") # aggregation en sommant la biomasse par pixel

        ras_ll <- project(ras_agg, "EPSG:4326") # latlon for leaflet & attention, reprojection raster modifie les valeurs dans les pixels
        values(ras_ll)[!is.na(values(ras_ll))] <- 1 # raster ne presente aucune valeur =0, sont représentes par des NA
        rast_tot <- c(rast_tot, ras_ll)
    }

    print(paste0("----------> ", species[i, 1], " DONE !"))
}

rast_stack <- rast(rast_tot)
rs <- sum(rast_stack, na.rm = TRUE) # calcul rs

x11()
terra::plot(rs, main = "Boulanger species RS")

writeRaster(rs, paste0("/home/claire/BDQC-GEOBON/data/EBVs/rs_Boulanger.tif"), overwrite = TRUE)

cmd <- paste0("s5cmd cp -acl public-read /home/claire/BDQC-GEOBON/data/EBVs/rs_Boulanger.tif s3://bq-io/acer/ebv/")
system(cmd)

#### Visualisation  ####
# -------------------- #

library(terra)
library(leaflet)
library(sf)
?aggregate
rast <- rast("/vsicurl/https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/rs_Boulanger.tif")


pal3 <- colorNumeric(c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"), values(rast), na.color = "transparent")

pal2 <- colorNumeric(rev(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")), values(rast), na.color = "transparent")

leaflet() %>%
    addTiles() %>%
    addRasterImage(rast)
