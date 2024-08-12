# Obtention de tendance temporelle de richesse specifiques par division spatiales du qc

library(terra)
library(sf)
library(mapview)
library(RCurl)
library(rmapshaper)
library(dplyr)
library(gdalcubes)
library(rstac)
library(stringr)

# CERQ avec 5 nibveaux
# Ecozones
ecoz <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecozones.gpkg")
ecoz <- ms_simplify(ecoz)
bbox <- st_bbox(ecoz)

# Ecoregions
ecor <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecoregions.gpkg")
ecor <- ms_simplify(ecor)

# Ecoprovinces
# Ecodistricts
ecod <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg")
ecod <- ms_simplify(ecod)


species <- read.table("INLA_species_list_final.txt", h = F)

#### Rstac & gdalcube ###
# Objectif: dev un workflow qui utilise directement les cartes de distribution dans le stac catalogue, calcul la richesse specifique et envoie le resultat, a nouveau, dans le stac catalogue


#### gdalcubes

# Connexion to the ACER stac catalog
catal <- stac("https://acer.biodiversite-quebec.ca/stac/")

## List collections dans ACER
collections <- catal %>%
    collections() %>%
    get_request() # get the list of collections

# Show the collection "oiseaux-nicheurs-qc"
obj <- catal %>%
    collections(collection_id = "oiseaux-nicheurs-qc") %>%
    items(limit = 30000) %>%
    get_request()

# Filtrer en fn des proprietes des items et creer une collection
range_coll <- gdalcubes::stac_image_collection(obj$features,
    asset_names = c("data"),
    property_filter = function(f) {
        f[["map"]] %in% c("range")
    }
) # ici f correspond a obj[['features']][[1]]$properties

range_vireo_soli_2017 <- stac_image_collection(obj$features,
    asset_names = c("data"),
    property_filter = function(f) {
        f[["map"]] %in% c("range") & f[["year"]] == 2017 & str_detect(f[["species"]], "Vireo solitarius") == TRUE
    }
)

# creation d'un cube
# To create a regular raster data cube from the image collection, we define the geometry of our target cube as a data cube view, using the cube_view() function. We define a simple overview, covering the full spatiotemporal extent of the imagery at 1km x 1km pixel size where one data cube cell represents a duration of one year. The provided resampling and aggregation methods are used to spatially reproject, crop, and rescale individual images and combine pixel values from many images within one year respectively. The raster_cube() function returns a proxy object, i.e., it returns immediately without doing any expensive computations.

cv <- cube_view(
    extent = extent(range_vireo_soli_2017),
    srs = "EPSG:6624",
    dt = "P1Y",
    dx = 10000,
    dy = 10000,
    aggregation = "median",
    resampling = "bilinear"
)
# Jumeler la collection et le cube_view pour créer un raster cube.
ras_cub <- raster_cube(range_vireo_soli_2017, cv)
plot(ras_cub, zlim = c(0, 10000)) # raster entirement noir




library(RColorBrewer)
raster_cube(range_vireo_soli_2017, cv) |>
    select_bands(c("data")) |>
    apply_pixel(function(v) {
        sum(v)
    }, names = "Species_richness") |>
    plot(zlim = c(0, 1), nbreaks = 10, col = brewer.pal(9, "YlGn"), key.pos = 1)



# EXAMPLE - Utilisez le jeu de données CHELSA sur les climatologies et créez une carte des moyennes pour les mois de juin, juillet et août 2010 à 2019
s_obj <- stac("https://io.biodiversite-quebec.ca/stac/")
bbox <- st_bbox(c(xmin = -76, xmax = -70, ymax = 54, ymin = 50), crs = st_crs(4326))
it_obj <- s_obj |>
    stac_search(collections = "chelsa-monthly", datetime = "2010-06-01T00:00:00Z/2019-08-01T00:00:00Z") |>
    get_request() |>
    items_fetch()

v <- cube_view(
    srs = "EPSG:32198", extent = list(
        t0 = "2010-06-01", t1 = "2019-08-31",
        left = bbox$xmin, right = bbox$xmax, top = bbox$ymax, bottom = bbox$ymin
    ),
    dx = 1000, dy = 1000, dt = "P10Y",
    aggregation = "mean",
    resampling = "bilinear"
)

for (i in 1:length(it_obj$features)) {
    names(it_obj$features[[i]]$assets) <- "data"
    it_obj$features[[i]]$assets$data$roles <- "data"
}
anames <- unlist(lapply(it_obj$features, function(f) {
    f["id"]
}))

st <- stac_image_collection(it_obj$features, asset_names = "data", property_filter = function(f) {
    f[["variable"]] == "tas" & (f[["month"]] %in% c(6, 7, 8))
})
c_cube <- raster_cube(st, v)
c_cube |> plot(col = heat.colors)

c_cube |> write_tif("~/", prefix = "chelsa-monthly", creation_options = list("COMPRESS" = "DEFLATE"))







# Sauvegarder le fichier résultant sur votre ordinateur

# lc_cube |> write_tif('~/',prefix='lc2',creation_options=list('COMPRESS'='DEFLATE'))

pal <- colorRampPalette(c("black", "darkblue", "red", "yellow", "white"))
lc_cube |> plot(zlim = c(0, 100), col = pal(10))


#### Sandbox ####
# summ_rs <- data.frame(matrix(ncol = 4, nrow = 0))
# x <- c("year", "poly_name", "poly_ID", "spe_rich")
# colnames(summ_rs) <- x

# for (year in 1992:2017) {
#     print(paste0("---> année en cours ", year))
#     path_list <- c()
#     # spe_rich <- 0
#     for (i in 1:length(species[, 1])) {
#         path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/", species[i, 1], "_range_", year, "_cog.tif")
#         print(path)
#         if (url.exists(path)) {
#             print("Path exists")
#             path_list <- c(path_list, paste0("/vsicurl_streaming/", path))
#         }
#     }
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






#     for (j in 56:nrow(ecod)) { # 1
#         cp <- terra::crop(stak, ecod[j, ])
#         rs_poly <- sum(colSums(cp[, -1], na.rm = T) > 0)
#         rw <- c(year, "ecod", j, rs_poly)
#         summ_rs <- rbind(summ_rs, rw)
#         print(paste0("poly # ", j, " done on ", nrow(ecod), " polygons"))
#     }
# } # WARNING, pb avec poly 55
# # }

# # mapview(ecoz[7, ])
# ras1 <- rast("/home/local/USHERBROOKE/juhc3201/Desktop/zonotrichia_leucophrys_range_2015.tif")
# ras2 <- rast("/home/local/USHERBROOKE/juhc3201/Desktop/zonotrichia_leucophrys_range_2015_COG.tif")

# x11()
# par(mfrow = c(1, 2))
# terra::plot(ras1)
# terra::plot(ras2)


# Sys.time()
# stak <- rast(path_list)
# Sys.time()
