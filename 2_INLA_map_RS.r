#### Prod carte RS from INLA maps (birds only)
# ----------------------------------------- #
library(terra)
library(gdalcubes)
library(RCurl) # for check if url exists - url.exists()

species <- read.csv("INLA_species_list_final.txt", header = F, col.names = c("species"))

# id_feat_Vince <- paste0(species_select, "_range_2017")
# paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/", id_feat_Vince, ".tif")

# rast1 <- rast("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/setophaga_cerulea_range_2017.tif")
# rast2 <- rast("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/regulus_satrapa_range_2017.tif")


for (year in 1993:2017) {
    print(paste0("--------------------> YEAR ", year))
    rast_tot <- c()
    for (i in 1:length(species$species)) {
        path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/oiseaux-nicheurs-qc/", species[i, 1], "_range_", year, ".tif")

        if (url.exists(path)) {
            map <- rast(path)
            map_ll <- project(map, "EPSG:4326") # lat lon for leaflet & ATTENTION, reproj modifie les valeurs dans les pixels
            rast_tot <- c(rast_tot, map_ll)
        }

        print(paste0("----------> ", species[i, 1], " DONE !"))
    }

    rast_stack <- rast(rast_tot)
    rs <- sum(rast_stack) # calcul rs
    names(rs) <- year # for keeping year information
    x11()
    terra::plot(rs, main = year)

    writeRaster(rs, paste0("/home/claire/BDQC-GEOBON/data/EBVs/rs_inla/rs_inla_", year, ".tif"), overwrite = TRUE)

    cmd <- paste0("s5cmd cp -acl public-read /home/claire/BDQC-GEOBON/data/EBVs/rs_inla/rs_inla_", year, ".tif s3://bq-io/acer/ebv/rs_inla/")
    system(cmd)
}

#### Visualisation rs_inla ####
# -------------------------- #

library(terra)
library(leaflet)
library(sf)

rast <- rast("/vsicurl/https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/rs_inla/rs_inla_1992.tif")


pal <- colorNumeric(c("#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506"), values(rast), na.color = "transparent")
pal3 <- colorNumeric(c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"), values(rast), na.color = "transparent")

pal2 <- colorNumeric(rev(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")), values(rast), na.color = "transparent")

leaflet() %>%
    addTiles() %>%
    addRasterImage(rast, colors = pal3, opacity = 0.8) %>%
    addLegend(pal = pal3, values = values(rast), title = "Richesse spécifique")
