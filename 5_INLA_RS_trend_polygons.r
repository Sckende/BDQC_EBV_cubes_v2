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


species <- read.table("INLA_species_list_final.txt", h = F)

summ_rs <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("year", "poly_name", "poly_ID", "spe_rich")
colnames(summ_rs) <- x

for (year in 1992:2017) {
    print(paste0("---> année en cours ", year))
    path_list <- c()
    spe_rich <- 0
    for (i in 1:length(species[, 1])) {
        path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/", species[i, 1], "_range_", year, "_cog.tif")
        print(path)
        if (url.exists(path)) {
            print("Path exists")
            path_list <- c(path_list, paste0("/vsicurl_streaming/", path))
        }
    }
    ### FAUT UTILISER GDALCUBE !!!!!
    synth_poly <- list()
    for (link in 1:length(path_list)) {
        stak <- rast(path_list[link])

        cp <- terra::extract(stak, ecoz)
        names(cp)[2] <- "occ"
        synth <- cp %>%
            group_by(ID) %>%
            summarize(max = max(occ, na.rm = T))

        synth$year <- year

        synth_poly[[link]] <- synth

        print(paste0("Year ", species[i, 1], " is done !"))






        # for (j in 1:nrow(ecoz)) {
        #     cp <- terra::crop(stak, ecoz[j, ])
        #     rs_poly <- sum(colSums(cp[, -1], na.rm = T) > 0)
        #     rw <- c(year, "ecoz", j, rs_poly)
        #     summ_rs <- rbind(summ_rs, rw)
        #     print(paste0("poly # ", j, " done"))
        # }
    }
}

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
df <- data.frame(
    id = character(),
    title = character(),
    description = character()
)
for (c in collections[["collections"]]) {
    df <- rbind(
        df,
        data.frame(
            id = c$id,
            title = c$title,
            description = c$description
        )
    )
}
df # table with collection names and their description

# Show the collection "oiseaux-nicheurs-qc"
obj <- catal %>%
    stac_search(collections = "oiseaux-nicheurs-qc") %>%
    get_request()


# obj %>% items_matched()
obj %>% items_length()

## See item properties of first item
obj[["features"]][[1]]$properties

## DL on local one item
# dl_ass <- obj %>% assets_download(assets_name = "zonotrichia_leucophrys_range_2017")


## Get one item and send it to STARS
library(stars) # good complement of TERRA packages for manipulating raster
lc1 <- read_stars(paste0("/vsicurl/", obj[["features"]][[2]]$assets$data$href),
    proxy = TRUE
) # by putting "/vsicurl/" and proxy = T it's just using the web, we avoid to read the raster on the computer memory => really faster !
plot(lc1)


input <- "acanthis_flammea"

id_feat <- paste(input, "_range_2017", sep = "")


acan <- stac("https://acer.biodiversite-quebec.ca/stac/") %>%
    collections("oiseaux-nicheurs-qc") %>%
    items(feature_id = id_feat) %>%
    get_request()

tif <- acan$assets$data$href

lc1 <- stars::read_stars(paste0("/vsicurl/", tif),
    proxy = TRUE
)
plot(lc1, axes = T)
