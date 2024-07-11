# Objectif: Pour chaque pixel de 1x1 inclut dans le Qc, obtenir la liste des espèces et la richesse spécifique
library(sf)
library(terra)
library(viridis)
library(dplyr)
library(stringr)
library(duckdb)
library(duckdbfs)
library(sfarrow)

data_path <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2"

# Local Atlas parquet

atlas_local <- function(parquet_file,
                        tblname = "atlas") {
    requireNamespace("duckdbfs")
    atlas <- duckdbfs::open_dataset(parquet_file, tblname = tblname)
    atlas
}

atlas <- atlas_local(
    paste0(data_path, "/atlas_2024-05-29.parquet"),
    "atlas"
)
dim(atlas)

# Check if obs only within Qc
atlas |>
    group_by(within_quebec) |>
    summarize(cnt = count())

# Conversion in sf object
# WARNING - VOIR si possible d'éviter cette conversion en travaillant directement avec un format geoparquet, qui permettrait des requêtes spatiales directement sur le fichier
colnames(atlas)

atlas2 <- atlas |> mutate_at(c("longitude", "latitude"), as.numeric)
atlas_sf <- atlas2 |>
    mutate(geometry = ST_Point(longitude, latitude)) |>
    to_sf(crs = 4326)

atlas_sf






# Conversion pix Qc from utm to lonlat
pix <- st_read("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale.gpkg")
st_crs(pix)
pixll <- st_transform(pix, crs = st_crs(4326))
head(pixll)
st_write(pixll, "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale_latlon.gpkg")




















# Association pixels-ecozones
# -------------------------- #
qc_ecozones <- st_read(paste0(data_path, "/qc_polygons/qc_ecozones.gpkg"))
names(qc_ecozones)
unique(qc_ecozones$ZONE_ID)
plot(st_geometry(qc_ecozones), col = viridis(length(unique(qc_ecozones$ZONE_ID))))
table(qc_ecozones$ECOZONE)


for (i in 1:length(st_geometry(qc_ecozones))) {
    poly <- qc_ecozones[i, ]
    print(poly$ZONE_NAME)
    wkt <- st_as_text(st_geometry(poly))

    pix <- st_read(paste0(data_path, "/qc_polygons/qc_grid_1x1km_finale.gpkg"),
        wkt_filter = wkt
    )

    pix$TYPE <- "ECOZONE"
    pix$ID_TYPE <- poly$ECOZONE
    pix$ECOZONE_NAME <- poly$ZONE_NAME

    if (i == 1) {
        st_write(
            pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecoz.gpkg")
        )
    } else {
        st_write(pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecoz.gpkg"),
            append = T
        )
    }
}

# d <- st_read(paste0(data_path, "/pix_infos/qc_pix_info_ecoz.gpkg"))

# Association pixels-ecoregions
# ------------------------------ #

qc_ecoregions <- st_read(paste0(data_path, "/qc_polygons/qc_ecoregions.gpkg"))
names(qc_ecoregions)
unique(qc_ecoregions$ECOREGION)
plot(st_geometry(qc_ecoregions), col = viridis(length(unique(qc_ecoregions$ECOREGION))))
plot(st_geometry(qc_ecoregions), border = "grey", add = T)
table(qc_ecoregions$ECOREGION)


for (i in 1:length(st_geometry(qc_ecoregions))) {
    poly <- qc_ecoregions[i, ]
    print(poly$REGION_NAM)
    wkt <- st_as_text(st_geometry(poly))

    pix <- st_read(paste0(data_path, "/qc_polygons/qc_grid_1x1km_finale.gpkg"),
        wkt_filter = wkt
    )

    pix$TYPE <- "ECOREGION"
    pix$ID_TYPE <- poly$ECOREGION
    pix$ECOZONE_NAME <- poly$REGION_NAM

    if (i == 1) {
        st_write(
            pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecor.gpkg")
        )
    } else {
        st_write(pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecor.gpkg"),
            append = T
        )
    }
}

# Association pixels-ecoprovinces
# ----------------------------- #

qc_ecoprovinces <- st_read(paste0(data_path, "/qc_polygons/qc_ecoprovinces.gpkg"))
names(qc_ecoprovinces)
unique(qc_ecoprovinces$ECOPROVINC)
plot(st_geometry(qc_ecoprovinces), col = viridis(length(unique(qc_ecoprovinces$ECOPROVINC))))
plot(st_geometry(qc_ecoprovinces), border = "grey", add = T)
table(qc_ecoprovinces$ECOPROVINC)


for (i in 1:length(st_geometry(qc_ecoprovinces))) {
    poly <- qc_ecoprovinces[i, ]
    print(poly$ECOPROVINC)
    wkt <- st_as_text(st_geometry(poly))

    pix <- st_read(paste0(data_path, "/qc_polygons/qc_grid_1x1km_finale.gpkg"),
        wkt_filter = wkt
    )

    pix$TYPE <- "ECOPROVINCE"
    pix$ID_TYPE <- poly$ECOPROVINC
    pix$ECOZONE_NAME <- poly$ECOPROVINC

    if (i == 1) {
        st_write(
            pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecop.gpkg")
        )
    } else {
        st_write(pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecop.gpkg"),
            append = T
        )
    }
}

# Association pixels-ecodistricts
# ------------------------------ #

qc_ecodistricts <- st_read(paste0(data_path, "/qc_polygons/qc_ecodistricts.gpkg"))
names(qc_ecodistricts)
unique(qc_ecodistricts$ECODISTRIC)
plot(st_geometry(qc_ecodistricts), col = viridis(length(unique(qc_ecodistricts$ECODISTRIC))))
plot(st_geometry(qc_ecodistricts), border = "grey", add = T)
table(qc_ecodistricts$ECODISTRIC)


for (i in 1:length(st_geometry(qc_ecodistricts))) {
    poly <- qc_ecodistricts[i, ]
    print(poly$ECODISTRIC)
    wkt <- st_as_text(st_geometry(poly))

    pix <- st_read(paste0(data_path, "/qc_polygons/qc_grid_1x1km_finale.gpkg"),
        wkt_filter = wkt
    )

    pix$TYPE <- "ECODISTRICT"
    pix$ID_TYPE <- poly$ECODISTRIC
    pix$ECOZONE_NAME <- poly$ECODISTRIC

    if (i == 1) {
        st_write(
            pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecod.gpkg")
        )
    } else {
        st_write(pix,
            paste0(data_path, "/pix_infos/qc_pix_info_ecod.gpkg"),
            append = T
        )
    }
}

## Traitement pour pixels dupliqués

pix_ecoz <- st_read(paste0(data_path, "/pix_infos/qc_pix_info_ecoz.gpkg"))
pix_ecop <- st_read(paste0(data_path, "/pix_infos/qc_pix_info_ecop.gpkg"))
pix_ecor <- st_read(paste0(data_path, "/pix_infos/qc_pix_info_ecor.gpkg"))
pix_ecod <- st_read(paste0(data_path, "/pix_infos/qc_pix_info_ecod.gpkg"))

length(st_geometry(pix_ecoz))
length(st_geometry(pix_ecop))
length(st_geometry(pix_ecor))
length(st_geometry(pix_ecod))



ecod_131 <- pix_ecod[pix_ecod$ID_TYPE == 131, ]
plot(st_geometry(ecod_131))
dup <- pix_ecod[duplicated(pix_ecod$ID), ]
plot(st_geometry(dup))






# dis <- qc_ecodistricts[, c(5, 6, 8)]
# names(dis)[1:2] <- c("ECODISTRICT_ID", "ECOREGION_ID")

# dis$ECOZONE_ID <- NA
# dis$ECOPROVINCE_ID <- NA

# for(i in 1:length(st_geometry(qc_ecozones))){
#     poly <- qc_ecozones[i, ]

#     sub <- dis[poly, ]

#     dis$ECOZONE_ID[dis$ECODISTRICT_ID %in% sub$ECODISTRICT_ID] <- poly$ECOZONE
# }

# for(i in 1:length(st_geometry(qc_ecoprovinces))){
#     poly <- qc_ecoprovinces[i, ]

#     sub <- dis[poly, ]

#     dis$ECOPROVINCE_ID[dis$ECODISTRICT_ID %in% sub$ECODISTRICT_ID] <- poly$ECOZONE
# }

#### Autre possibilite : rasterizer toutes les occurrences par année par espèce puis extraire par wkt

## phaese de test
atlas_local <- function(parquet_file,
                        tblname = "atlas") {
    requireNamespace("duckdbfs")
    atlas <- duckdbfs::open_dataset(parquet_file, tblname = tblname)
    atlas
}

atlas <- atlas_local("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_2024-05-29.parquet") |>
    filter(within_quebec == "t")

atlas_sp <- atlas |>
    filter(group_fr == "Mammifères") |>
    group_by(valid_scientific_name) |>
    summarize(cnt = n()) |>
    collect()

alces <- atlas |>
    filter(valid_scientific_name == "Alces americanus") |>
    group_by(year_obs) |>
    summarize(cnt = n()) |>
    select()

alces_year <- alces |>
    group_by(year_obs) |>
    summarize(cnt = n())

# From tibble to sf obj
alces_df <- as.data.frame(alces)
alces_sf <- st_as_sf(alces_df,
    coords = c("longitude", "latitude"),
    crs = 4326
)
alces
unique(alces_sf$year_obs)

alces_sf_2000 <- alces_sf[alces_sf$year_obs == "2000", ]
plot(st_geometry(alces_sf_2000))
summary(as.numeric(alces_df$latitude))
dim(alces_df[alces_df$year_obs == "2000" & as.numeric(alces_df$latitude) == -90, ])


#### Prod carte RS from INLA maps (birds only)
# ----------------------------------------- #
library(terra)
library(gdalcubes)
library(RCurl) # for check if url exists - url.exists()

species <- read.csv("INLA_species_list_final.txt", header = F, col.names = c("species"))

id_feat_Vince <- paste0(species_select, "_range_2017")
paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/", id_feat_Vince, ".tif")

rast1 <- rast("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/setophaga_cerulea_range_2017.tif")
rast2 <- rast("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/regulus_satrapa_range_2017.tif")


for (year in 1992:2017) {
    print(paste0("--------------------> YEAR ", year))
    rast_tot <- c()
    for (i in 1:length(species$species)) {
        path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/", species[i, 1], "_range_", year, ".tif")

        if (url.exists(path)) {
            rast_tot <- c(rast_tot, rast(path))
        }

        print(paste0("----------> ", species[i, 1], " DONE !"))
    }

    rast_stack <- rast(rast_tot)
    rs <- sum(rast_stack) # calcul rs
    names(rs) <- year # for keeping year information
    x11()
    terra::plot(rs, main = year)

    writeRaster(rs, paste0("/home/claire/BDQC-GEOBON/data/EBVs/rs_inla/rs_inla_", year, ".tif"), overwrite = TRUE)

    cmd <- paste0("s5cmd cp -acl public-read /home/claire/BDQC-GEOBON/data/EBVs/rs_inla/rs_inla_", year, ".tif s3://bq-io/acer/ebv/rs_inla/")
    system(cmd)
}

#### Visualisation rs_inla ####
# -------------------------- #

library(terra)
library(leaflet)

rast <- rast("/home/claire/BDQC-GEOBON/data/EBVs/rs_inla/rs_inla_1992.tif")
crs(rast)
rast2 <- project(rast, "EPSG:4326") # lat lon for leaflet

pal <- colorNumeric(c("#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506"), values(rast2), na.color = "transparent")
pal3 <- colorNumeric(c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"), values(rast2), na.color = "transparent")

pal2 <- colorNumeric(rev(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")), values(rast2), na.color = "transparent")

leaflet() %>%
    addTiles() %>%
    addRasterImage(rast2, colors = pal3, opacity = 0.8) %>%
    addLegend(pal = pal3, values = values(rast2), title = "Richesse spécifique")
