library(gdalcubes)
library(rstac)
library(knitr)
library(stringr)
library(stars)
library(terra)
library(dplyr)

# Connexion au catalogue STAC ACER
acer <- stac("https://acer.biodiversite-quebec.ca/stac/")
io <- stac("https://io.biodiversite-quebec.ca/stac/")

# Liste des collections contenues dans le catalogue
collections <- acer |>
    collections() |>
    get_request()

# Obtenir le nom des collections et leur description
df <- data.frame(
    id = character(),
    title = character(),
    description = character()
)
for (c in collections[["collections"]]) {
    df <- rbind(df, data.frame(id = c$id, title = c$title, description = c$description))
}
kable(df)

# Recuperer une collection specifique, ici "oiseaux-nicheurs-qc"
obj <- acer %>%
    collections(collection_id = "oiseaux-nicheurs-qc") %>%
    items(limit = 30000) %>%
    get_request()

# Voir les proprietes du premier item (ou couche)
obj[["features"]][[1]]$properties

# Filtrer les items d'une collection a partir de leurs proprietes
stac_items <- acer |>
    ext_filter(
        collection %in% c("oiseaux-nicheurs-qc") &&
            `map` %in% c("range") &&
            `year` == 1992
    ) |>
    get_request() |>
    items_fetch() # item_fetch() pour recuperer tous les items au travers des differentes paginations

# Proprietes des items
stac_items$features
# lien vers un des items
stac_items[["features"]][[100]]$assets$data$href

# Liste des noms des items dans la collection
for (i in 1:length(stac_items$features)) {
    print(stac_items$features[[i]]$id)
    # print(stac_items$features[[i]]$links)
}

# visualisation de l'item #100 avec le package STARS
it <- read_stars(paste0("/vsicurl/", stac_items[["features"]][[80]]$assets$data$href), proxy = TRUE)
plot(it)

# Calcul de la richesse specifique a l'echelle du qc
href_ls <- c()
for (i in 1:length(stac_items$features)) {
    href <- paste0("/vsicurl/", stac_items$features[[i]]$assets$data$href)
    href_ls <- c(href_ls, href)
}

ras <- rast(href_ls)
ras2 <- sum(ras)
plot(ras2)

#### FUNCTION pour extraire richesse specifique dans le temps a l'echelle du QC ####
sr_trend <- function(catalog = "acer", collection = "oiseaux-nicheurs-qc", years = NULL) {
    stac_cat <- stac(paste0("https://", catalog, ".biodiversite-quebec.ca/stac/"))
    # spe_num_df <- data.frame(year = numeric(), spe_rich = numeric())
    spe_num_df <- tibble(year = numeric(), spe_rich = numeric(), spe_ls = list())


    if (is.null(years)) {
        print("Computation for full time range - from 1992 to 2015 - quiet long ")
        years <- 1992:2017
    }

    for (i in years) {
        print(i)
        stac_items <- stac_cat |>
            ext_filter(
                collection == collection &&
                    `map` %in% c("range") &&
                    `year` == {{ i }}
            ) |>
            get_request() |>
            items_fetch()

        spe_num <- stac_items |>
            items_length()

        spe_ls <- NULL
        for (j in 1:spe_num) {
            spe <- stac_items$features[[j]]$properties$species
            spe_ls <- c(spe_ls, spe)
        }


        spe_num_df <- rbind(spe_num_df, tibble(year = i, spe_rich = spe_num, spe_ls = list(spe_ls)))
    }

    spe_num_df
}


sr_trend()
sr_trend(years = 1995)
sr_trend(years = 2011:2015)


#### fonction pour calcul de RS dans le temps pour un polygone donne ####
# see https://rdrr.io/cran/gdalcubes/man/extract_geom.html
ecoz <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecozones.gpkg")
ecoz <- rmapshaper::ms_simplify(ecoz)

ecod <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg")
ecod <- rmapshaper::ms_simplify(ecod)
poly <- ecod[10:11, ]

plot(st_geometry(ecod))
plot(st_geometry(poly), col = "grey", add = T)

# vireos <- rast("/home/claire/desktop/vireo_data_test/vireo_single_file.tif")
range <- rast("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/inla_multiband/2017_range_multiband.tif")
names(range)

plot(range)


# ext <- extract(range2017, vect(poly))
# max_ext <- apply(ext, 2, max)
# names(max_ext)[max_ext > 0][-1]

sr_spat_temp_trend <- function(catalog = "acer", collection = "oiseaux-nicheurs-qc", years = NULL, polygon = NULL) {
    if (!is(polygon, "sf") & !is.null(polygon)) {
        stop("Polygon must be a sf object")
    }
    polygon$poly_id <- 1:nrow(polygon)

    if (is.null(years)) {
        print("Computation for full time range - from 1992 to 2015 - quiet long ")
        years <- 1992:2017
    }
    if (min(years) < 1992 | max(years) > 2017) {
        stop("Years have to be from 1992 to 2017")
    }
    # Empty df for storing the result
    spe_num_df <- tibble(year = numeric(), spe_rich = numeric(), spe_ls = list(), poly_id = numeric())

    for (p in 1:nrow(polygon)) {
        poly <- polygon[p, ]
        print(paste0("Computation for polygon  number ", p, " on ", nrow(polygon)))

        for (y in years) {
            # print(y)
            # Retrieve the multi-band raster !!!! A MODIFIER QUAND FILES DANS ACER STAC !!!
            # multi <- rast(paste0("/home/claire/desktop/data/multiband/", y, "_range_multiband_renamed.tif")) # yoga
            multi <- rast(paste0("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/inla_multiband/", y, "_range_multiband_renamed.tif")) # xps


            # extract raster values in polygon for each band
            ext <- extract(multi, vect(poly))
            max_ext <- apply(ext, 2, max)
            spe_ls <- names(max_ext)[max_ext > 0][-1]
            spe_num <- length(spe_ls)

            # result
            # spe_num_df <- rbind(spe_num_df, tibble(year = y, spe_rich = spe_num, spe_ls = list(spe_ls), poly_id = p, geom = polygon[p, "geom"]))
            spe_num_df <- add_row(spe_num_df, year = y, spe_rich = spe_num, spe_ls = list(spe_ls), poly_id = p)
        }
    }

    final <- left_join(spe_num_df, polygon[, c("poly_id", "geom")], by = c("poly_id"))
    final
}

# Zone de test #
sr_spat_temp_trend(polygon = ecod)
sr_spat_temp_trend(years = 2000, polygon = ecod)
rs <- sr_spat_temp_trend(years = c(1992, 2003, 2017), polygon = ecod)
rs |> print(n = 303)

rs_ls <- split(rs, rs$poly_id)
sr_spat_temp_trend(years = 2019, polygon = poly)

ecoz_rs <- sr_spat_temp_trend(polygon = ecoz)
ecoz_rs_ls <- split(ecoz_rs, ecoz_rs$poly_id)

x11()
par(mfrow = c(3, 3))
lapply(ecoz_rs_ls, function(x) {
    plot(x$year, x$spe_rich, type = "l")
})


#### Fonction de Guillaume ###

library(gdalcubes)
library(rstac)
library(sf)
gdalcubes_options(parallel = TRUE)

ecod <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg")
ecod <- rmapshaper::ms_simplify(ecod)
ecod6622 <- st_transform(ecod, "EPSG:6622")
poly <- ecod6622[70, ]
bb <- st_bbox(poly)

s_obj <- stac("https://acer.biodiversite-quebec.ca/stac/")
it_obj <- s_obj |>
    stac_search(collections = "oiseaux-nicheurs-qc", limit = 1000) |>
    post_request() |>
    items_fetch()

v <- cube_view(
    extent = list(
        left = bb["xmin"], right = bb["xmax"],
        bottom = bb["ymin"], top = bb["ymax"], t0 = "1992-01-01", t1 = "1992-01-01"
    ),
    # extent = list(
    #     left = -500000, right = -200000,
    #     bottom = 250000, top = 400000, t0 = "1992-01-01", t1 = "1992-01-01"
    # ),
    srs = "EPSG:6622", dx = 10000, dy = 10000, dt = "P1Y", aggregation = "max"
)

# Test de cube view avec un extend de polygone
# poly <- ecoz[5, ]
# poly_bbox <- st_bbox(poly)

# v <- cube_view(
#     extent = list(
#         left = poly_bbox[1], right = poly_bbox[3],
#         bottom = poly_bbox[2], top = poly_bbox[4], t0 = "1992-01-01", t1 = "1992-01-01"
#     ),
#     srs = "EPSG:6622", dx = 10000, dy = 10000, dt = "P1Y", aggregation = "max"
# )

# v2 <- cube_view(srs = "EPSG:6622", dx = 10000, dy = 10000, dt = "P1Y", aggregation = "max", resampling = "average", extent = list(left = bbox["xmin"] - 1000, right = bbox["xmax"] + 1000,
#         bottom = bbox["ymin"] - 1000, top = bbox["ymax"] + 1000, t0 = "1992-01-01", t1 = "1992-01-01"))

# gdalcubes_options(parallel = 16)
# cb <- raster_cube(tc, v2) |>
#     filter_geom(poly$geom)
# plot(cb)

n <- data.frame(species = character(), max = numeric())
j <- 0
for (i in 1:items_length(it_obj)) {
    tmp_it <- it_obj
    tmp_it$features <- it_obj$features[[i]]
    if (tmp_it$features$properties$map == "pocc" & tmp_it$features$properties$year == "1992") {
        j <- j + 1
        tc <- stac_image_collection(it_obj$features[i], asset_names = c("data"))
        cube <- raster_cube(tc, v)
        cube |> write_tif("/home/local/USHERBROOKE/juhc3201/Downloads/GIS/tmp")
        arr <- cube |>
            apply_pixel("data>0.25") |>
            reduce_space("max(band1)") |> # plot possible avec valeurs par chaque pas de temps (si il y en a de defini)
            as_array()
        n[j, "species"] <- tmp_it$features$properties$species
        n[j, "max"] <- arr[1]
    }
}

















# https://appelmar.github.io/ogh2021/tutorial.html

library(rstac)
library(gdalcubes)
library(stars)
library(tmap)
library(sf)

# area of interest
ecoz <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecozones.gpkg")
ecoz <- rmapshaper::ms_simplify(ecoz)
poly <- ecoz[6, ]

tmap_mode("view")
tm_shape(st_geometry(poly)) + tm_polygons()

st_crs(poly)

# We aim at generating a cloud-free composite image of our study area for June, 2018 and we use the rstac package to find suitable Sentinel-2 images. However, to use the bbox argument of the corresponding function stac_search() for spatial filtering, we first need to derive and transform the bounding box to latitude / longitude (WGS84) values, for which we use the st_bbox() and st_transform() functions.

bbox <- st_bbox(poly)
bbox

st_as_sfc(bbox) |>
    st_transform("EPSG:4326") |>
    st_bbox() -> bbox_wgs84
bbox_wgs84

s <- stac("https://acer.biodiversite-quebec.ca/stac/")

items <- s |>
    stac_search(
        collections = "oiseaux-nicheurs-qc",
        bbox = c(
            bbox_wgs84["xmin"], bbox_wgs84["ymin"],
            bbox_wgs84["xmax"], bbox_wgs84["ymax"]
        ),
        datetime = "2000-01-01/2000-12-31",
        limit = 500
    ) |>
    post_request()
items
