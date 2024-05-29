# Objectif: obtenir differentes echelles spatiales au sein du Quebec: ecozone, ecoregion, ecoprovince, ecodistrict, les différents niveaux de cadres ecologiques (niv 01 a 05) et un grille de pixels de 1x1 km applicable pour le calcul de la richesse specifique a partir des donnees brutes de Atlas
# Libraries
library(sf)
library(mapview)

# canadian ecozones, ecoprovinces, ecoregions, ecodistricts
# data downloaded from https://sis.agr.gc.ca/siscan/nsdb/ecostrat/gis_data.html
local_path <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_regions"
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
local_path_save <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons"
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

qc_grid <- st_intersection(grid, qc) # super long, test sur Narval: /home/ccjuhasz/projects/def-dgravel/ccjuhasz/QC_in_a_CUBE/version_2
