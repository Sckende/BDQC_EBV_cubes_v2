library(stringr)
library(terra)
year <- Sys.getenv("YEARBLOC") # retrieve the year variable from bash script
print(year)

file_n <- strsplit(list.files("/home/claire/desktop/data/", pattern = ".tif"), "_")

name_ls <- lapply(file_n, function(x) {
    paste(x[1], x[2], sep = "_")
})

rs <- rast(paste0("/home/claire/desktop/data/multiband/", year, "_range_multiband.tif"))

names(rs) <- unlist(name_ls)

writeRaster(rs, paste0("/home/claire/desktop/data/multiband/", year, "_range_multiband_renamed.tif"), overwrite = TRUE)
