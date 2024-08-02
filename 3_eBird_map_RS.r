#### Prod carte RS from eBird maps (birds only)
# ----------------------------------------- #
# la liste d'especes correspond aux oiseaux nicheurs du Qc
# identique a celle utilisee par V. Bellavance
# une seule carte de RS produite car pas de tendance annuelle pour eBird

library(terra)
library(gdalcubes)
library(RCurl) # for check if url exists - url.exists()

species <- read.csv("INLA_species_list_final.txt", header = F, col.names = c("species"))

rast_tot <- c()
for (i in 1:length(species$species)) {
    path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/eBird_maps/", species[i, 1], "_range.tif")

    if (url.exists(path)) {
        ras <- rast(path)
        ras_ll <- project(ras, "EPSG:4326") # latlon for leaflet & attention, reprojection raster modifie les valeurs dans les pixels
        values(ras_ll)[values(ras_ll) > 0 & !is.na(values(ras_ll))] <- 1
        rast_tot <- c(rast_tot, ras_ll)
    }

    print(paste0("----------> ", species[i, 1], " DONE !"))
}

rast_stack <- rast(rast_tot)
rs <- sum(rast_stack) # calcul rs

x11()
terra::plot(rs, main = "eBird RS")

writeRaster(rs, paste0("/home/claire/BDQC-GEOBON/data/EBVs/rs_ebird.tif"), overwrite = TRUE)

cmd <- paste0("s5cmd cp -acl public-read /home/claire/BDQC-GEOBON/data/EBVs/rs_ebird.tif s3://bq-io/acer/ebv/")
system(cmd)

#### Visualisation  ####
# -------------------- #

library(terra)
library(leaflet)
library(sf)

rast <- rast("/vsicurl/https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/rs_ebird.tif")

pal3 <- colorNumeric(c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"), values(rast), na.color = "transparent")

pal2 <- colorNumeric(rev(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")), values(rast), na.color = "transparent")

leaflet() %>%
    addTiles() %>%
    addRasterImage(rast, colors = pal3, opacity = 0.8) %>%
    addLegend(pal = pal3, values = values(rast), title = "Richesse spécifique")

#### Pour envoi dans stac catalogue ####
# conversion en Quebec Albers (projection conique equal area)

r <- rast("/vsicurl/https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/rs_ebird.tif")
repj <- terra::project(r, "EPSG:6623")

writeRaster(repj, "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/io-data/bqio/data/rs_ebird.tif", overwrite = T)
