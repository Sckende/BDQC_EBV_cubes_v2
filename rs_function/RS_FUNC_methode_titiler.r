library(jsonlite)
library(geojsonsf)
library(httr2)
library(rstac)
library(sf)
library(terra)
library(rmapshaper)
library(dplyr)

# Qc unique polygon
qc <- ms_simplify(st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg"))
qc1 <- st_as_sf(as.data.frame(st_combine(qc)))
qc4326 <- st_transform(qc1, "EPSG:4326")
st_geometry(qc4326) <- "geom"

poly <- qc[10, ]
# poly4326 <- st_transform(poly, "EPSG:4326")
# geoj_sf <- sf_geojson(poly4326)



sr_extract <- function(catalog = "acer", collection = "oiseaux-nicheurs-qc", years = NULL, polygon = NULL) {
    # catalog
    cata <- stac(paste0("https://", catalog, ".biodiversite-quebec.ca/stac/"))

    # Empty df for storing the result
    spe_num_df <- tibble(year = numeric(), spe_rich = numeric(), spe_ls = list())

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

    for (p in 1:nrow(polygon)) {
        poly <- polygon[p, ]
        print(paste0("Computation for polygon  number ", p, " on ", nrow(polygon)))

        poly_geoj <- sf_geojson(poly)

        for (y in years) {
            it_obj <- cata |>
                ext_filter(
                    collection == {{ collection }} &&
                        `map` %in% c("range") &&
                        `year` == {{ y }}
                ) |>
                get_request() |>
                items_fetch()



            # first step for extracting data
            n2 <- data.frame(species = character(), max = numeric())
            print(paste0("----------> ", items_length(it_obj), " species found for ", y))

            for (i in 1:items_length(it_obj)) {
                # print(it_obj$features[[i]]$assets$data$href)
                print(paste0(i, "/", items_length(it_obj)))
                req <- request("https://tiler.biodiversite-quebec.ca/cog/statistics") |>
                    # req_body_raw(geojson_sf(toJSON(it_obj$features[[i]])) |>
                    req_body_raw(poly_geoj) |> # WARNING ! Polygon has to be in geojson & in 4326 !!!
                    # st_set_crs(st_crs("EPSG:6622")) |>
                    # st_transform("EPSG:4326") |>
                    # sf_geojson()) |>
                    req_url_query(url = it_obj$features[[i]]$assets$data$href, categorical = "TRUE", max_size = 1024) |>
                    req_perform() |>
                    resp_body_json()

                n2[i, "species"] <- it_obj$features[[i]]$properties$species
                n2[i, "max"] <- req$features[[1]]$properties$statistics$b1$max
            }
            spe_rich <- sum(n2$max, na.rm = T)
            spe_ls <- sub(" ", "_", (tolower(n2$species[n2$max != 0])))

            spe_num_df <- add_row(spe_num_df, year = y, spe_rich = spe_rich, spe_ls = list(spe_ls))
        }
    }

    spe_num_df
}


#### test zone ####

# RS pour le Qc en 1992
test1b <- sr_extract(years = c(2000))
# RS pour le Qc pour toutes les annees (1992:2017)
test2 <- sr_extract()
# RS pour deux polygones en 2000
test3 <- sr_extract(polygon = qc[c(5, 30), ], years = 2000)
# RS pour deux polygones en 1998 & 2010
test4 <- sr_extract(polygon = qc[c(10, 70), ], years = c(1998, 2010))




spe <- sub(" ", "_", tolower(spe_ls))
nospe <- sub(" ", "_", tolower(n2$species[n2$max == 0]))
x11()
par(mfrow = c(4, 3))
for (i in 1:length(nospe)) {
    lk <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/oiseaux-nicheurs-qc/", nospe[i], "_range_2000.tif")
    terra::plot(rast(lk))
    plot(st_geometry(st_transform(poly4326, "EPSG:6622")), add = T)
}
