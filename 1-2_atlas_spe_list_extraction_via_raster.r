# Methode via raster directement pour extraction de la liste des especes de Atlas par pixel appartenant au Quebec
# code source original -> François Rousseu

library(duckdbfs)
library(dplyr)
library(sf)
library(geodata)
library(rmapshaper)
library(data.table)

epsg <- 6623

can <- gadm("CAN", level = 1, path = getwd()) |> st_as_sf()
qc <- can[can$NAME_1 %in% c("Québec"), ] |>
    st_transform(crs = epsg) |>
    ms_simplify(0.01)

# x11(); par(mfrow = c(1, 2))
# plot(ms_simplify(st_geometry(can), 0.01)); plot(st_geometry(qc))

# Chargement des fonctions necessaires pour utilisation d'Atlas
source("https://object-arbutus.cloud.computecanada.ca/bq-io/atlas/parquet/bq-atlas-parquet.R")

# Requete Atlas en remote
atlas <- atlas_remote(parquet_date = "2024-07-16") # version geoparquet
colnames(atlas)
atlas |>
    count(group_en) |>
    print(n = 100)

taxa <- atlas |>
    distinct(group_en) |>
    pull() # pull() permet de retrieve le result sous format de vecteur

# Creation du raster pour la creation des cartes de richesse specifique
res <- 10000
r <- rast(resolution = res, extent = st_bbox(st_buffer(qc, res)))

# Selection des occurrences dans Atlas et conversion dans une projection equal area (Quebec Albers)

for (tax in taxa) {
    occs <- atlas |> filter(group_en == tax)

    years <- occs |>
        distinct(year_obs) |>
        pull()

    for (y in years) {
        occs_y <- occs |>
            filter(year_obs == y) |>
            mutate(geom = ST_Point(as.numeric(longitude), as.numeric(latitude))) |>
            to_sf(crs = 4326) |>
            collect() |>
            st_transform(epsg)

        # Conserve les occs uniquement dans le Qc
        occs_y <- occs_y[qc, ]

        ids <- cells(r, vect(occs_y))

        nbsp <- cbind(occs_y, ids) |>
            st_drop_geometry() |>
            _[, c("valid_scientific_name", "cell")] |> # "_" represente le resultat precedent et est utilise pour faire de l'indexation
            as.data.table() |>
            unique() |> # inverse d'un duplicated, considere toutes les colonnes dans le datatable
            _[, .(n = .N), by = "cell"] # syntaxe de datatable - permet de faire somme des especes unique par cell

        r[as.integer(nbsp$cell)] <- nbsp$n

        x11()
        plot(r, main = y)
        plot(st_geometry(qc), add = TRUE)
    }
}
