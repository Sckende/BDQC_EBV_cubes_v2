# Aires de répartition des mammifères terrestres, des reptiles, des amphibiens et des poissons d'eau douce
#  sources : https://www.donneesquebec.ca/recherche/dataset/aires-de-repartition-faune


# amphibiens : 21 espèces
# reptiles : 17 espèces
# mammifères terrestres : 69 espèces
# poissons d'eau douce et migrateurs : 118 espèces

library(sf)
library(terra)

amph <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/data/Aires_repartition_amphibiens.gpkg")
amph_l <- split(amph, amph$NOM_FRANCA)
mapview::mapview(amph_l)
rept <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/data/Aires_repartition_reptiles.gpkg")
rept_l <- split(rept, rept$NOM_FRANCA)
mapview::mapView(rept_l)
mamm <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/data/Aires_repartition_MT.gpkg")
mamm_l <- split(mamm, mamm$NOM_FRANCA)
mapview::mapview(mamm_l)
poiss <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/data/Aires_repartition_poisson_eau_douce.gpkg")
poiss_l <- split(poiss, poiss$NOM_FRANCA)
mapview::mapview(poiss[113, ])

# Creation du raster de ref pour la rasterisation des polygones

# recup du bbox du Quebec en UTM
ecoz <- st_read("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/qc_polygons/sf_eco_poly/qc_ecozones.gpkg")
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
mapview::mapview(ref_ras)
terra::plot(ref_ras)

# rasterization test
# gdalcubes::extent(ref_ras) <- extent(poly)
pol_test <- mamm[1, ]
rp <- rasterize(pol_test, ref_ras, "NOM_FRANCA")

mapview::mapview(rp)

# rasterization looooooop
melccfp <- list(amph_l, rept_l, mamm_l, poiss_l)
lapply(melccfp, function(x) { #
    lapply(x, function(y) {
        rp <- rasterize(y, ref_ras, "NOM_FRANCA")
        print(y$NOM_FRANCA)

        # traitement pour commande bash pour envoyer raster dans bucket s3
        dfy <- as.data.frame(y)
        col_group <- stringr::str_which(names(dfy), "GRAND_GROU")
        group <- dfy[1, col_group]
        spe <- stringr::str_replace_all(unique(dfy$NOM_ANGLA), " ", "_")
        col_date <- stringr::str_which(names(dfy), "DATE")
        year <- dfy[1, col_date]

        # send the raster to s3
        path_s3 <- paste0("s3://bq-io/acer/melccfp/", group, "/melccfp_", spe, "_", year, ".tif")
        path <- paste0("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/rasterisation/melccfp_", spe, "_", year, ".tif")


        writeRaster(rp,
            path,
            filetype = "GTiff",
            gdal = "COMPRESS=DEFLATE"
        )
        bash <- paste("AWS_ACCESS_KEY_ID=NJBPPQZX7PFUBP1LH8B0 AWS_SECRET_ACCESS_KEY=DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd cp -acl public-read", path, path_s3, sep = " ")
        system(bash)
        system(paste("rm", path))
    })
})
