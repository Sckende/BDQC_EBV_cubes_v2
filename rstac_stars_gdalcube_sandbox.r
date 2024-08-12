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
ecod <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg")
ecod <- rmapshaper::ms_simplify(ecod)
poly <- ecod[70, ]

plot(st_geometry(ecod))
plot(st_geometry(poly), col = "grey", add = T)

# vireos <- rast("/home/claire/desktop/vireo_data_test/vireo_single_file.tif")
range2017 <- rast("/home/claire/desktop/vireo_data_test/2017_range_single_file.tif")


ext <- extract(range2017, vect(poly))
max_ext <- apply(ext, 2, max)
names(max_ext)[max_ext > 0][-1]




# def de la bbox
poly_bbox <- poly |>
    sf::st_transform(32198) |>
    sf::st_bbox()

poly_bbox_geojson <- rstac::cql2_bbox_as_geojson(poly_bbox)
poly_bbox_geojson

# filtrage des collection
poly_item <- acer |>
    ext_filter(
        collection %in% c("oiseaux-nicheurs-qc") &&
            `map` %in% c("range") &&
            `species` == "Zonotrichia leucophrys" &&
            `year` == 1992
    ) |>
    get_request() |>
    items_fetch()





# Pour optimisation et calcul selon des polygones
# ==> see stars package https://r-spatial.github.io/stars/
# ==> see gdalcube package https://github.com/appelmar/gdalcubes
# fusion des deux pour une solution !

f1 <- read_stars(href_ls[1])
f2 <- read_stars(href_ls[2])
new <- c(f1, f2)

r_stars <- read_stars(c(href_ls[1], href_ls[2], href_ls[3]), proxy = TRUE, along = "attr") # argument along is for determining how several arrays are combined
r_stars |>
    slice(index = 1, along = "attr") |>
    plot()

# multiband file with gdalcubes

vireo_col <- create_image_collection("/home/claire/desktop/vireo_data_test/vireo_single_file.tif")








vireos <- brick("/home/claire/desktop/vireo_data_test/vireo_single_file.tif")



names(vireos)
band








lc1

bbox <- st_bbox(c(xmin = -76, xmax = -70, ymax = 54, ymin = 50), crs = st_crs(4326))
lc2 <- lc1 |> st_crop(bbox)

pal <- colorRampPalette(c("black", "darkblue", "red", "yellow", "white"))
plot(lc2, breaks = seq(0, 100, 10), col = pal(10))

for (i in 1:length(it_obj$features)) {
    it_obj$features[[i]]$assets$data$roles <- "data"
}

st <- stac_image_collection(it_obj$features, asset_names = c("data"), property_filter = function(f) {
    f[["class"]] %in% c("1", "2", "3", "4")
}, srs = "EPSG:4326")
st


bbox <- st_bbox(c(xmin = -483695, xmax = -84643, ymin = 112704, ymax = 684311), crs = st_crs(32198))

v <- cube_view(srs = "EPSG:32198", extent = list(
    t0 = "2000-01-01", t1 = "2000-01-01",
    left = bbox$xmin, right = bbox$xmax, top = bbox$ymax, bottom = bbox$ymin
), dx = 1000, dy = 1000, dt = "P1D", aggregation = "sum", resampling = "mean")

lc_cube <- raster_cube(st, v)
lc_cube |> plot(zlim = c(0, 100), col = pal(10))


it_obj <- s_obj |>
    collections("accessibility_to_cities") |>
    items() |>
    get_request() |>
    items_fetch()
v <- cube_view(srs = "EPSG:32198", extent = list(
    t0 = "2015-01-01", t1 = "2015-01-01",
    left = bbox$xmin, right = bbox$xmax, top = bbox$ymax, bottom = bbox$ymin
), dx = 1000, dy = 1000, dt = "P1D", aggregation = "mean", resampling = "bilinear")
for (i in 1:length(it_obj$features)) {
    it_obj$features[[i]]$assets$data$roles <- "data"
}
st <- stac_image_collection(it_obj$features)
lc_cube <- raster_cube(st, v)
lc_cube |> plot(col = heat.colors)

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

it_obj[["features"]][[1]]$properties





create_image_collection()
