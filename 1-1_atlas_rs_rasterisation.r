library(sf)
library(terra)
library(arrow)
library(stringr)


# Creation du raster de reference
ecoz <- st_read("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_ecozones.gpkg")
bb <- st_bbox(ecoz)

ref_ras <- rast(
    xmin = bb[1],
    xmax = bb[3],
    ymin = bb[2],
    ymax = bb[4],
    resolution = 1000
)
crs(ref_ras) <- "epsg:32198"
values(ref_ras) <- 1
ref_ras_ll <- terra::project(ref_ras, "epsg:4326")

# Loop pour rasterisation de chaq parquet cree par groupe taxo par annee
pq_path <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/rasterisation/"
ras_path <- "/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/raster/"
pq_ls <- list.files(pq_path, full.names = TRUE)

for (i in 1:length(pq_ls)) {
    pq_path <- pq_ls[i]
    y <- read_parquet(pq_path)
    y$geometry <- st_as_sfc(y$geometry)
    y <- st_as_sf(y,
        sf_column_name = "geometry",
        crs = 4326
    )

    # rasterisation
    rp <- rasterize(y, ref_ras_ll, "spe_rich", touches = TRUE, fun = "mean") # for fun argument, mean value of RS if two polygons are touched by a pixel
    info <- str_sub(pq_path, (str_locate(pq_path, "//")[1, 2]) + 1, (str_locate(pq_path, "_rs.parquet")[1, 1]) - 1)
    print(paste0("-------------------------------------> ", info))

    names(rp) <- info

    write_ras_path <- paste0(ras_path, info, "_rs_atlas_2024-07-16.tif")
    path_s3 <- paste0("s3://bq-io/acer/ebv/rs_atlas/", info, "_rs_atlas_2024-07-16.tif")
    writeRaster(rp,
        write_ras_path,
        filetype = "GTiff",
        gdal = "COMPRESS=DEFLATE",
        overwrite = TRUE
    )
    bash <- paste("AWS_ACCESS_KEY_ID=xxx AWS_SECRET_ACCESS_KEY=xxx S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd cp -acl public-read", write_ras_path, path_s3, sep = " ")
    system(bash)
    system(paste("rm", write_ras_path))
}
