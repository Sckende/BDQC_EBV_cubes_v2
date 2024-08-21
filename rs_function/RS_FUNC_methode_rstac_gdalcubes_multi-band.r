# --------> TO DO: envoyer gtif multiband dans stac catalogue pour tester

library(gdalcubes)
library(rstac)
library(sf)
library(terra)
library(dplyr)
library(rmapshaper)
gdalcubes_options(parallel = TRUE)

#### fonction pour calcul de RS dans le temps pour un polygone donne avec COGgtif multiband ####
# ------------------------------------------------------------------------------------------- #
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

#### Fonction ####
# -------------- #

# Qc unique polygon
qc <- ms_simplify(st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg"))
qc1 <- st_as_sf(as.data.frame(st_combine(qc)))
qc32198 <- st_transform(qc1, "EPSG:32198")
qc4326 <- st_transform(qc, "EPSG:4326")
st_geometry(qc32198) <- "geom"

poly <- qc[10, ]

sr_spat_temp_trend <- function(catalog = "acer", collection = "oiseaux-nicheurs-qc", years = NULL, polygon = NULL) {
    # Polygon conditions
    if (is.null(polygon)) {
        polygon <- qc32198 # polygone du Qc par defaut
    } else if (!is(polygon, "sf")) {
        stop("Polygon must be a sf object")
    } else {
        polygon <- st_transform(polygon, "EPSG:32198") # !!!!! WARNING ! To change depending on the raster crs !
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
            multi <- rast(paste0("/home/claire/desktop/data/multiband/", y, "_range_multiband_renamed.tif")) # yoga
            # multi <- rast(paste0("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/inla_multiband/", y, "_range_multiband_renamed.tif")) # xps


            # extract raster values in polygon for each band
            ext <- extract(multi, vect(poly), touches = TRUE)
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
# RS pour le Qc en 1992
test1b <- sr_spat_temp_trend(years = c(1992))
# RS pour le Qc pour toutes les annees (1992:2017)
test2b <- sr_spat_temp_trend()
# RS pour deux polygones en 2000
test3b <- sr_spat_temp_trend(polygon = qc[c(5, 30), ], yers = 2000)
# RS pour deux polygones en 1998 & 2010
test4b <- sr_spat_temp_trend(polygon = qc[c(10, 70), ], yers = c(1998, 2010))

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
