library(sf)
library(viridis)
data_path <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons"
qc_ecozones <- st_read(paste0(data_path, "/qc_ecozones.gpkg"))
names(qc_ecozones)
unique(qc_ecozones$ZONE_ID)
plot(st_geometry(qc_ecozones), col = viridis(length(unique(qc_ecozones$ZONE_ID))))

for (i in length(st_geometry(qc_ecozones))) {
    poly <- qc_ecozones[i, ]
    wkt <- st_as_text(st_geometry(obj))

    pix <- st_read(paste0(data_path, "/qc_grid_1x1km_finale.gpkg"),
        wkt_filter = wkt
    )

    pix$ECOZONE <- poly$ECOZONE
    pix$ECOZONE_NAME <- poly$ZONE_NAME
}

plot(st_geometry(qc_ecozones)[qc_ecozones$ECOZONE == 7])
