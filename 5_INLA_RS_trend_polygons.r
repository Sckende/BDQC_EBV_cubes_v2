# Obtention de tendance temporelle de richesse specifiques par division spatiales du qc

library(terra)
library(sf)
library(mapview)
library(RCurl)
library(rmapshaper)
library(dplyr)
library(gdalcubes)
library(rstac)

# CERQ avec 5 nibveaux
# Ecozones
ecoz <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecozones.gpkg")
ecoz <- ms_simplify(ecoz)

# Ecoregions
ecor <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecoregions.gpkg")
ecor <- ms_simplify(ecor)

# Ecoprovinces
# Ecodistricts
ecod <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg")
ecod <- ms_simplify(ecod)


species <- read.table("INLA_species_list_final.txt", h = F)

summ_rs <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("year", "poly_name", "poly_ID", "spe_rich")
colnames(summ_rs) <- x

for (year in 1992:2017) {
    print(paste0("---> année en cours ", year))
    path_list <- c()
    # spe_rich <- 0
    for (i in 1:length(species[, 1])) {
        path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/", species[i, 1], "_range_", year, "_cog.tif")
        print(path)
        if (url.exists(path)) {
            print("Path exists")
            path_list <- c(path_list, paste0("/vsicurl_streaming/", path))
        }
    }
    # ### FAUT UTILISER GDALCUBE !!!!!
    # synth_poly <- list()
    # for (link in 1:length(path_list)) {
    #     stak <- rast(path_list[link])

    #     cp <- terra::extract(stak, ecoz)
    #     names(cp)[2] <- "occ"
    #     synth <- cp %>%
    #         group_by(ID) %>%
    #         summarize(max = max(occ, na.rm = T))

    #     synth$year <- year

    #     synth_poly[[link]] <- synth

    #     print(paste0("Species ", species[i, 1], " is done !"))






    for (j in 56:nrow(ecod)) { # 1
        cp <- terra::crop(stak, ecod[j, ])
        rs_poly <- sum(colSums(cp[, -1], na.rm = T) > 0)
        rw <- c(year, "ecod", j, rs_poly)
        summ_rs <- rbind(summ_rs, rw)
        print(paste0("poly # ", j, " done on ", nrow(ecod), " polygons"))
    }
} # WARNING, pb avec poly 55
# }

# mapview(ecoz[7, ])
ras1 <- rast("/home/local/USHERBROOKE/juhc3201/Desktop/zonotrichia_leucophrys_range_2015.tif")
ras2 <- rast("/home/local/USHERBROOKE/juhc3201/Desktop/zonotrichia_leucophrys_range_2015_COG.tif")

x11()
par(mfrow = c(1, 2))
terra::plot(ras1)
terra::plot(ras2)


Sys.time()
stak <- rast(path_list)
Sys.time()


#### gdalcubes

# Connexion to the stac catalog
catal <- stac("https://acer.biodiversite-quebec.ca/stac/")

## List collections
collections <- catal %>%
    collections() %>%
    get_request() # get the list of collections

## Show collections and descriptions
# library(knitr)
# df<-data.frame(id=character(),title=character(),description=character())
# for (c in collections[['collections']]){
#   df<-rbind(df,data.frame(id=c$id,title=c$title,description=c$description))
# }
# kable(df)

# Show the collection "oiseaux-nicheurs-qc"
obj <- catal %>%
    collections(collection_id = "oiseaux-nicheurs-qc") %>%
    items(limit = 30000) %>%
    get_request() # à remettre à jour quand on saura depasser la hard limit de l'api a 10000

# Filtrer en fn des proprietes des item et creer une collection
range_coll <- stac_image_collection(obj$features, asset_names = c("data"), property_filter = function(f) {
    f[["map"]] %in% c("range")
}) # ici f correspond a obj[['features']][[1]]$properties

range_coll <- stac_image_collection(obj$features, asset_names = c("data"), property_filter = function(f) {
    f[["map"]] %in% c("range") & f[["year"]] == 2017
})

# creation d'un cube
ecoz <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecozones.gpkg")

bbox <- st_bbox(ecoz[6, ]) # a mettre à jour avec les polygones ecoz
bbox <- st_bbox(ecoz)

v <- cube_view(
    srs = "EPSG:32198",
    extent = list(
        t0 = "2017-01-01",
        t1 = "2017-01-01",
        left = bbox$xmin,
        right = bbox$xmax,
        top = bbox$ymax,
        bottom = bbox$ymin
    ),
    dx = 1024,
    dy = 1024,
    dt = "P1Y",
    aggregation = "sum",
    resampling = "mean"
)

# Jumeler la collection et le cube_view pour créer un raster cube.

lc_cube <- raster_cube(range_coll, v)
dim(lc_cube)
print(lc_cube)
summary(as.data.frame(lc_cube))

# Sauvegarder le fichier résultant sur votre ordinateur

# lc_cube |> write_tif('~/',prefix='lc2',creation_options=list('COMPRESS'='DEFLATE'))

pal <- colorRampPalette(c("black", "darkblue", "red", "yellow", "white"))
lc_cube |> plot(zlim = c(0, 100), col = pal(10))
