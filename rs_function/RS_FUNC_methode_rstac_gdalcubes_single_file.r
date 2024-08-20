# ---------------------------> TO DO: construire la fonction
library(gdalcubes)
library(rstac)
library(sf)
library(rmapshaper)
library(dplyr)

gdalcubes_options(parallel = TRUE)

# Qc unique polygon
qc <- ms_simplify(st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg"))
qc1 <- st_as_sf(as.data.frame(st_combine(qc)))
qc6622 <- st_transform(qc1, "EPSG:6622")
st_geometry(qc6622) <- "geom"

sr_spat_temp_trend <- function(catalog = "acer", collection = "oiseaux-nicheurs-qc", years = NULL, polygon = NULL) {
    # Polygon conditions
    if (is.null(polygon)) {
        polygon <- qc6622 # polygone du Qc par defaut
    } else if (!is(polygon, "sf")) {
        stop("Polygon must be a sf object")
    } else {
        polygon <- st_transform(polygon, "EPSG:6622")
    }
    polygon$poly_id <- 1:nrow(polygon)

    # Year conditions
    if (is.null(years)) {
        print("Computation for full time range - from 1992 to 2015 - quiet long ")
        years <- 1992:2017
    }
    if (min(years) < 1992 | max(years) > 2017) {
        stop("Years have to be from 1992 to 2017")
    }

    # Catalog
    clg <- stac(paste0("https://", catalog, ".biodiversite-quebec.ca/stac/"))

    # Empty df for storing the result
    spe_num_df <- tibble(year = numeric(), spe_rich = numeric(), spe_ls = list(), poly_id = numeric())

    for (y in years) {
        # Items depending on year
        it_obj <- clg |>
            ext_filter(
                collection == {{ collection }} &&
                    `map` %in% c("range") &&
                    `year` == {{ y }}
            ) |>
            get_request() |>
            items_fetch()

        for (p in 1:nrow(polygon)) {
            print(paste0("Computation for polygon  number ", p, " on ", nrow(polygon), " - Year ", y))

            # get the bbox of the polygon
            bb <- st_bbox(polygon[p, ])
            # define the cube dimension
            v <- cube_view(
                extent = list(
                    left = bb["xmin"], right = bb["xmax"],
                    bottom = bb["ymin"], top = bb["ymax"],
                    t0 = paste0(y, "-01-01"), t1 = paste0(y, "-01-01")
                ),
                srs = "EPSG:6622",
                dx = 10000, dy = 10000,
                dt = "P1Y",
                aggregation = "max"
            )

            # dataframe to retrieve data
            n <- data.frame(species = character(), max = numeric())
            for (i in 1:items_length(it_obj)) {
                print(paste0("-----> ", i, " species on ", items_length(it_obj)))
                # for (i in 1:20) { # !!!!! A CHANGER POUR LE BENCHMARK DES FONCTIONS
                # tmp_it <- it_obj

                # tmp_it$features <- it_obj$features[[i]]

                tc <- create_image_collection(
                    c(paste0("/vsicurl/", it_obj$features[[i]]$assets$data$href)),
                    date_time = c(paste0(y, "-01-01")),
                    band_names = "data"
                )
                cube <- raster_cube(tc, v)
                arr <- cube |>
                    filter_geom(polygon[p, ]$geom) |>
                    reduce_space("max(data)") |>
                    as_array()

                n[i, "species"] <- it_obj$features[[i]]$properties$species
                n[i, "max"] <- arr[1]
            }
            # info summary
            spe_num <- sum(n$max, na.rm = T)
            spe_ls <- n$species[n$max != 0]
            spe_num_df <- add_row(spe_num_df, year = y, spe_rich = spe_num, spe_ls = list(spe_ls), poly_id = p)
        }
    }
    final <- left_join(spe_num_df, polygon[, c("poly_id", "geometry")], by = c("poly_id"))
    final
}

# zone de test #
sr_spat_temp_trend(years = 2000)
