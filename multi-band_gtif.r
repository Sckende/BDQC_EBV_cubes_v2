library(stringr)
library(terra)
year <- Sys.getenv("YEARBLOC") # retrieve the year variable from bash script
print(year)

file_n <- strsplit(list.files("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/dl_acer/", pattern = ".tif"), "_")

name_ls <- lapply(file_n, function(x) {
    paste(x[1], x[2], sep = "_")
})

rs <- rast(paste0("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/inla_multiband/", year, "_range_multiband.tif"))

names(rs) <- unlist(name_ls)

writeRaster(rs, paste0("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/inla_multiband/", year, "_range_multiband_renamed.tif"), overwrite = TRUE)
