# Packages
library(jsonlite)
library(geojsonsf)
library(httr2)
library(rstac)
library(sf)
library(terra)
library(rmapshaper)
library(dplyr)

# Data necessaire pour faire tourner la fonction et pour les exemples
# Qc unique polygon
qc <- ms_simplify(st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecodistricts.gpkg"))
qc1 <- st_as_sf(as.data.frame(st_combine(qc)))
qc4326 <- st_transform(qc1, "EPSG:4326")
st_geometry(qc4326) <- "geom"


#' @param catalog nom du catalogue stac a utiliser - par defaut acer
#' @param collection nom de la collection utilisee pour extraire la richesse specifique - par defaut oiseaux-nicheurs-qc
#' @param years vecteur contenant les annees d'interet - si null, l'extraction est faite pour le range temporel maximal, cad 1992 à 2017
#' @param polygon objet sf contenant un ou plusieurs polygones d'interet pour l'extraction de la richesse specifique - si null, l'extraction est faite à l'échelle du Québec au complet
#' @return une valeur de richesse specifique par année par pollygone
#' @examples 
#' RS pour le Qc en 2000
#'test1 <- sr_extract(years = c(2000))
#' 
#' RS pour le Qc pour toutes les annees (1992:2017)
#'test2 <- sr_extract()
#' 
#' RS pour deux polygon pour toutes les annes
#'test3 <- sr_extract(polygon = qc[c(30, 90), ])
#' visualisation de test3
#' par(mfrow = c(1, 2))
#' plot(st_geometry(qc)); plot(st_geometry(unique(test3$geom)), add = T, col = c("darkorange", "seagreen"))
#' plot(test3$year, test3$spe_rich, xlab = "annee", ylab = "richesse specifique", type = "b", col = "darkorange")
#' 
#' RS pour deux polygones en 2000
#'test4 <- sr_extract(polygon = qc[c(5, 30), ], years = 2000)
#' 
#' RS pour deux polygones en 1998 & 2010
#'test5 <- sr_extract(polygon = qc[c(10, 70), ], years = c(1998, 2010))

#' @export
sr_extract <- function(catalog = "acer", collection = "oiseaux-nicheurs-qc", type = c("spe_rich", "area"), years = NULL, polygon = NULL, species = NULL) {
    # catalog
    cata <- stac(paste0("https://", catalog, ".biodiversite-quebec.ca/stac/"))
    
    # year condition
    if (is.null(years)) {
        years <- 1992:2017
    }

    # Empty df for storing the result
    spe_num_df <- tibble(year = numeric(), spe_rich = numeric(), spe_ls = list(), polygon_id = numeric())

    # SELECTION type = spe_rich
    if(type == "spe_rich"){
            # Polygon conditions
    if (is.null(polygon)) {
        polygon <- qc4326 # polygone du Qc par defaut
    } else if (!is(polygon, "sf")) {
        stop("Polygon must be a sf object")
    } else {
        polygon <- st_transform(polygon, "EPSG:4326") # devra etre modifiable en fn du crs des cartes utilisees
    }
    polygon$polygon_id <- 1:nrow(polygon)

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
                print(paste0(i, "/", items_length(it_obj)))
                req <- request("https://tiler.biodiversite-quebec.ca/cog/statistics") |>
                    req_body_raw(poly_geoj) |>
                    req_url_query(url = it_obj$features[[i]]$assets$data$href, categorical = "TRUE", max_size = 1024) |>
                    req_perform() |>
                    resp_body_json()

                n2[i, "species"] <- it_obj$features[[i]]$properties$species
                n2[i, "max"] <- req$features[[1]]$properties$statistics$b1$max
            }
            spe_rich <- sum(n2$max, na.rm = T)
            spe_ls <- sub(" ", "_", (tolower(n2$species[n2$max != 0])))

            spe_num_df <- add_row(spe_num_df, year = y, spe_rich = spe_rich, spe_ls = list(spe_ls), polygon_id = p)
        }
    }

    final <- left_join(spe_num_df, polygon[, c("polygon_id", "geom")], by = c("polygon_id"))
    final
    } else {
        if(is.null(species)){
        stop("Species must be indicated")
        }

    spe_area_df <-  tibble(year = numeric(), species = character(), spe_area_km2 = numeric())

    for (spe in species){
        print(spe)
        for (y in years) {
            it_obj <- cata |>
                ext_filter(
                    collection == {{ collection }} &&
                    `map` %in% c("range") &&
                    `year` == {{ y }} &&
                    `species`== {{ spe }}
                ) |>
                get_request() |>
                items_fetch()
        
            if(length(it_obj$feature) == 0) {
                spe_area = 0
            } else {
                map <- rast(it_obj$features[[1]]$assets$data$href)
                map_pj <- project(map, "EPSG:6623", method = "near") # conversion en proj equal area

                spe_area <- sum(values(map_pj)) * res(map_pj)[1] * res(map_pj)[1] / 1000000 # calcul de l'aire a partir de la res du raster en km2
            }

            spe_area_df <- add_row(spe_area_df, year = y, species = spe, spe_area_km2 = spe_area)
        }
    }
        
 

    final <- spe_area_df
    }

final

}

# ----- #
# Obtenir la liste des especes pour lesquelles des cartes sont contenues dans ACER
# ----- #

# --> catalogue
cata <- stac(paste0("https://", catalog, ".biodiversite-quebec.ca/stac/"))
# --> liste des items
it_obj <- cata |>
        ext_filter(
            collection == {{ "oiseaux-nicheurs-qc" }} &&
            `map` %in% c("range")
            ) |>
            get_request() |>
            items_fetch()
# --> obtenir la liste des especes            
spe <- c()
for(i in 1:500){ # 500 est random number exageremment grand pour etre sur d'avoir toutes les especes
    spe <- c(spe, it_obj[['features']][[i]]$properties$species)
}
spe_ls <- unique(spe)
spe_ls

# ----- #
# tests
# ----- #

test1 <- sr_extract(catalog = "acer", collection = "oiseaux-nicheurs-qc", type = "area", years = 1998, species = "Falcipennis canadensis")
test2 <- sr_extract(catalog = "acer", collection = "oiseaux-nicheurs-qc", type = "area", years = c(1998, 2005), species = "Falcipennis canadensis")
test3 <- sr_extract(catalog = "acer", collection = "oiseaux-nicheurs-qc", type = "area", years = NULL, species = "Falcipennis canadensis")
test4 <- sr_extract(catalog = "acer", collection = "oiseaux-nicheurs-qc", type = "area", years = NULL, species = "Catharus bicknelli")
test5 <- sr_extract(catalog = "acer", collection = "oiseaux-nicheurs-qc", type = "area", years = NULL, species = c("Aquila chrysaetos", "Catharus bicknelli", "Tyrannus tyrannus"))
test6 <- sr_extract(catalog = "acer", collection = "oiseaux-nicheurs-qc", type = "area", years = NULL, species = spe_ls[1:12])



plot(test4$year, test4$spe_area, type = "b")
par(mfrow = c(2, 2))
lapply(split(test5, test5$species), function(x){
    plot(x = x$year, y = x$spe_area_km2, main = unique(x$species), type = "b")
})
par(mfrow = c(3, 4))
lapply(split(test6, test6$species), function(x){
    plot(x = x$year, y = x$spe_area, main = unique(x$species), type = "b")
})
#### TO DO 
# 1 - calcul reelle de la superficie en fonction de la proj
map_pj <- project(map, "EPSG:6623", method = "near")
sum(values(map_pj)) * 10020.97 * 10020.97 / 1000000

sum(values(map)) * 10020.97 * 10020.97 / 1000000
