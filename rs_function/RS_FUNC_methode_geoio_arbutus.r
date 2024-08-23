library(jsonlite) #*
library(geojsonsf) #*
library(httr2) #*
library(rstac) #*
library(sf) #*
# library(terra)
library(rmapshaper) #*
library(dplyr) #*

# Qc unique polygon
qc <- ms_simplify(st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg"))
qc1 <- st_as_sf(as.data.frame(st_combine(qc)))
qc4326 <- st_transform(qc1, "EPSG:4326")
st_geometry(qc4326) <- "geom"



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
        polygon <- qc4326 # polygone du Qc par defaut
    } else if (!is(polygon, "sf")) {
        stop("Polygon must be a sf object")
    } else {
        polygon <- st_transform(polygon, "EPSG:4326") # !!!!! WARNING ! To change depending on the raster crs !
    }
    polygon$poly_id <- 1:nrow(polygon)
    for (y in years) {
        print(paste0("-----> year ", y))
        it_obj <- cata |>
            ext_filter(
                collection == {{ collection }} &&
                    `map` %in% c("range") &&
                    `year` == {{ y }}
            ) |>
            get_request() |>
            items_fetch()

        for (p in 1:nrow(polygon)) {
            polyg <- polygon[p, ]
            print(paste0("Computation for polygon  number ", p, " on ", nrow(polygon)))

            body <- list()
            body$cog_urls <- data.frame(ref = character(), url = character())
            # geoj <- geojson_sf(toJSON(it_obj$features[[1]], auto_unbox = TRUE)) |>
            #     st_set_crs(st_crs("EPSG:6622")) |>
            #     st_transform("EPSG:4326")

            geoj <- polyg |>
                st_transform("EPSG:4326")
            geoj$place <- "Area"
            body$geojson <- geoj |> sf_geojson(simplify = FALSE)

            # bboxsf <- geojson_sf(toJSON(it_obj$features[[1]])) |>
            #     st_set_crs(st_crs("EPSG:6622")) |>
            #     st_transform("EPSG:4326") |>
            #     sf_geojson()
            for (i in 1:items_length(it_obj)) {
                # tmp_it$features <- it_obj$features[[i]]
                # if (tmp_it$features$properties$map == "occ" & tmp_it$features$properties$year == "1992") {
                body$cog_urls <- rbind(body$cog_urls, data.frame(ref = it_obj$features[[i]]$properties$species, url = it_obj$features[[i]]$assets$data$href))
                # }
            }
            js <- toJSON(body, auto_unbox = TRUE)
            js <- gsub('\\\\"', '"', js)
            js <- gsub('\\\\"', '"', js)
            js <- gsub('\\\\"', '"', js)
            js <- gsub('\"geojson\":\"', '\"geojson\":', js)
            js <- gsub(']}\"}', "]}}", js)
            req <- request("https://geoio.biodiversite-quebec.ca/geojson_stats_many_urls/") |>
                req_body_raw(js) |>
                req_perform() |>
                resp_body_json()

            n3 <- data.frame(species = names(req), max = sapply(req, FUN = function(a) {
                a$Area_1$b1$max
            }))
            spe_rich <- sum(n3$max, na.rm = T)
            spe_ls <- sub(" ", "_", (tolower(n3$species[n3$max != 0])))

            spe_num_df <- add_row(spe_num_df, year = y, spe_rich = spe_rich, spe_ls = list(spe_ls), poly_id = p)
        }
    }
    spe_num_df
}


#### test zone ####
# RS pour le Qc en 1992
test1 <- sr_extract(years = c(1992))
# RS pour le Qc en 1992 et 1 polygone
test1_1 <- sr_extract(years = c(1992), polygon = qc[70, ])
# RS pour le Qc pour toutes les annees (1992:2017)
test2 <- sr_extract()
test2 |> print(n = 100)
# RS pour deux polygones en 2000
test3 <- sr_extract(polygon = qc[c(5, 30), ], years = 2000)
# RS pour deux polygones en 1998 & 2010
test4 <- sr_extract(polygon = qc[c(10, 70), ], years = c(1998, 2010))
