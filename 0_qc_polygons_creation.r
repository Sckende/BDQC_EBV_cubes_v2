# Objectif: obtenir differentes echelles spatiales au sein du Quebec: ecozone, ecoregion, ecoprovince, ecodistrict, les différents niveaux de cadres ecologiques (niv 01 a 05) et un grille de pixels de 1x1 km applicable pour le calcul de la richesse specifique a partir des donnees brutes de Atlas
# Libraries
library(sf)
library(mapview)

# canadian ecozones, ecoprovinces, ecoregions, ecodistricts
# data downloaded from https://sis.agr.gc.ca/siscan/nsdb/ecostrat/gis_data.html
# local_path <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_regions" # XPS
local_path <- "/home/claire/BDQC-GEOBON/data/QUEBEC_regions" # YOGA

ecoZo <- st_read(paste0(local_path, "/Ecozones")) # 25 polygones
ecoPro <- st_read(paste0(local_path, "/Ecoprovinces")) # 68 polygones
ecoReg <- st_read(paste0(local_path, "/Ecoregions")) # 218 polygones
ecoDis <- st_read(paste0(local_path, "/Ecodistricts")) # 1025 polygones
mapview(list(ecoZo, ecoPro, ecoReg, ecoDis))

# Quebec map
# data informations https://www.donneesquebec.ca/recherche/fr/dataset/cadre-ecologique-de-reference
# zip file accessible to https://stqc380donopppdtce01.blob.core.windows.net/donnees-ouvertes/Cadre_ecologique_reference/CERQ_SHP.zip
# missing step here is the conversion from shapefile to geopackage file
qc <- st_union(st_read(paste0(local_path, "/sf_CERQ_SHP/QUEBEC_CR_NIV_01.gpkg")))

# Crop canadian shapefile with Qc
can_spat <- list(ecoZo, ecoPro, ecoReg, ecoDis)
qc_eco <- lapply(can_spat, function(x) {
    x_conv <- st_transform(x, crs = st_crs(qc))
    inter <- st_intersection(x_conv, qc)
})

# Local save
# local_path_save <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons" # XPS
local_path_save <- paste0(local_path, "/sf_eco_poly") # XPS

n <- c("ecozones", "ecoprovinces", "ecoregions", "ecodistricts")
for (i in 1:4) {
    st_write(
        qc_eco[[i]],
        paste0(local_path_save, "/qc_", n[i], ".gpkg")
    )
}

# 1x1 km grid creation
cell_size <- c(1000, 1000) # en m, soit 1kmx1km

grid <- st_make_grid(qc,
    cellsize = cell_size,
    square = TRUE
)
# st_write(
#     grid,
#     paste0(local_path_save, "/qc_grid_1x1km_initiale.gpkg")
# )

# retrieve the wkt for qc
qc_wkt <- st_as_text(st_geometry(qc))

qc_1x1km <- st_read(paste0(local_path_save, "/qc_grid_1x1km_initiale.gpkg"),
    wkt_filter = qc_wkt
)

qc_1x1km$ID <- 1:length(st_geometry(qc_1x1km))

# x11(); plot(st_geometry(qc_1x1km)) # laborieux
st_write(
    qc_1x1km,
    paste0(local_path_save, "/qc_grid_1x1km_finale.gpkg")
)
