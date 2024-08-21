library(rstac)
library(terra)
library(sf)
library(dplyr)
library(rmapshaper)

# ecodistrict
qc <- ms_simplify(st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg"))

# Qc unique polygon
qc1 <- st_as_sf(as.data.frame(st_combine(qc)))
qc32198 <- st_transform(qc1, "EPSG:32198")


sr_extract <- function(catalog = "acer", collection = "oiseaux-nicheurs-qc", years = NULL, polygon = NULL) {
    # catalog
    cata <- stac(paste0("https://", catalog, ".biodiversite-quebec.ca/stac/"))

    # Empty df for storing the result
    spe_num_df <- tibble(year = numeric(), spe_rich = numeric(), spe_ls = list(), poly_id = numeric())

    # year condition
    if (is.null(years)) {
        years <- 1992:2017
    }

    # Polygon conditions
    if (is.null(polygon)) {
        polygon <- qc32198 # polygone du Qc par defaut
    } else if (!is(polygon, "sf")) {
        stop("Polygon must be a sf object")
    } else {
        polygon <- st_transform(polygon, "EPSG:32198") # !!!!! WARNING ! To change depending on the raster crs !
    }
    polygon$poly_id <- 1:nrow(polygon)


    for (y in years) {
        it_obj <- cata |>
            ext_filter(
                collection == {{ collection }} &&
                    `map` %in% c("range") &&
                    `year` == {{ y }}
            ) |>
            get_request() |>
            items_fetch()

        href <- vector()
        spe <- vector()
        for (l in 1:items_length(it_obj)) {
            href <- c(href, it_obj$features[[l]]$assets$data$href)
            spe <- c(spe, sub(" ", "_", tolower(it_obj$features[[l]]$properties$species)))
        }

        ras <- rast(href)
        names(ras) <- spe
        for (p in 1:nrow(polygon)) {
            poly <- polygon[p, ]
            print(paste0("Computation for polygon  number ", p, " on ", nrow(polygon)))
            # extract raster values in polygon for each band
            ext <- extract(ras, vect(poly), touches = TRUE)
            max_ext <- apply(ext, 2, max)
            spe_ls <- names(max_ext)[max_ext > 0][-1]
            spe_rich <- length(spe_ls)

            # result
            spe_num_df <- add_row(spe_num_df, year = y, spe_rich = spe_rich, spe_ls = list(spe_ls), poly_id = p)
        }
    }
    spe_num_df
}

#### test zone ####
# RS pour le Qc en 1992
test1 <- sr_extract(years = c(1992))
# RS pour le Qc pour toutes les annees (1992:2017)
test2 <- sr_extract()
# RS pour deux polygones en 2000
test3 <- sr_extract(polygon = qc[c(5, 30), ], years = 2000)
# RS pour deux polygones en 1998 & 2010
test4 <- sr_extract(polygon = qc[c(10, 70), ], years = c(1998, 2010))
