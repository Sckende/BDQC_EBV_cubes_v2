# Obtention de tendance temporelle de richesse specifiques par division spatiales du qc

library(terra)
library(sf)
library(mapview)

list.files("/home/claire/BDQC-GEOBON/data/QUEBEC_regions/sf_CERQ_SHP/")
# CERQ avec 5 nibveaux
# Ecozones
# Ecoregions
# Ecoprovinces
# Ecodistricts

# uniquement pour Yoga
local_path <- "/home/claire/BDQC-GEOBON/data/QUEBEC_regions/"
ecoZo <- st_read(paste0(local_path, "/Ecozones")) # 25 polygones
ecoPro <- st_read(paste0(local_path, "/Ecoprovinces")) # 68 polygones
ecoReg <- st_read(paste0(local_path, "/Ecoregions")) # 218 polygones
ecoDis <- st_read(paste0(local_path, "/Ecodistricts")) # 1025 polygones
mapview(list(ecoZo, ecoPro, ecoReg, ecoDis))

qc <- st_union(st_read(paste0(local_path, "/sf_CERQ_SHP/QUEBEC_CR_NIV_01.gpkg")))
qc_ll <- st_transform(qc, st_crs(ecoZo))

# Crop canadian shapefile with Qc
can_spat <- list(ecoZo, ecoPro, ecoReg, ecoDis)
qc_eco <- lapply(can_spat, function(x) {
    inter <- st_intersection(x, qc_ll)
})
