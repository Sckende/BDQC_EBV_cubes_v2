# ----------> TO DO: integrer filtrage par polygones + construire la fonction

library(jsonlite)
library(geojsonsf)
library(httr2)
library(rstac)
library(sf)
library(terra)
library(rmapshaper)

# Qc unique polygon
qc <- ms_simplify(st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg"))
qc1 <- st_as_sf(as.data.frame(st_combine(qc)))
qc4326 <- st_transform(qc1, "EPSG:4326")
st_geometry(qc32198) <- "geom"

poly <- qc[10, ]
poly4326 <- st_transform(poly, "EPSG:4326")

# catalog
acer <- stac("https://acer.biodiversite-quebec.ca/stac/")

it_obj <- acer |>
    ext_filter(
        collection == collection &&
            `map` %in% c("range") &&
            `year` == {{ 2000 }}
    ) |>
    get_request() |>
    items_fetch()

n2 <- data.frame(species = character(), max = numeric())

for (i in 1:items_length(it_obj)) {
    # tmp_it <- it_obj
    # tmp_it$features <- it_obj$features[[i]]

    print(it_obj$features[[i]]$assets$data$href)
    req <- request("https://tiler.biodiversite-quebec.ca/cog/statistics") |>
        req_body_raw(geojson_sf(toJSON(it_obj$features[[i]])) |>
            st_set_crs(st_crs("EPSG:6622")) |>
            st_transform("EPSG:4326") |>
            crop(vect(poly))
            sf_geojson()) |>
        req_url_query(url = it_obj$features[[i]]$assets$data$href, categorical = "TRUE") |>
        req_perform() |>
        resp_body_json()

    n2[i, "species"] <- tmp_it$features$properties$species
    n2[i, "max"] <- req$features[[1]]$properties$statistics$b1$max
}
n2
sum(n2$max)
