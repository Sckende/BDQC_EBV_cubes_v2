# Aires de répartition des mammifères terrestres, des reptiles, des amphibiens et des poissons d'eau douce
#  sources : https://www.donneesquebec.ca/recherche/dataset/aires-de-repartition-faune


# amphibiens : 21 espèces
# reptiles : 17 espèces
# mammifères terrestres : 69 espèces
# poissons d'eau douce et migrateurs : 118 espèces

#### Packages
library(sf)
library(terra)

#### Custom functions
accentless <- function(s) {
    chartr(
        "áéóūáéíóúÁÉÍÓÚýÝàèìòùÀÈÌÒÙâêîôûÂÊÎÔÛãõÃÕñÑäëïöüÄËÏÖÜÿçÇ",
        "aeouaeiouAEIOUyYaeiouAEIOUaeiouAEIOUaoAOnNaeiouAEIOUycC",
        s
    )
}

#### Data
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

#### Creation du raster de ref pour la rasterisation des polygones
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
# mapview::mapview(ref_ras)
# terra::plot(ref_ras)

#### rasterization looooooop
melccfp <- list(amph_l, rept_l, mamm_l, poiss_l)
lapply(melccfp, function(x) { #
    lapply(x, function(y) {
        # traitement pour commande bash pour envoyer raster dans bucket s3
        dfy <- as.data.frame(y)
        col_group <- stringr::str_which(names(dfy), "GRAND_GROU")
        group <- dfy[1, col_group]
        group <- tolower(accentless(group))
        spe <- stringr::str_replace_all(unique(dfy$NOM_FRANCA), " ", "-")
        spe <- tolower(accentless(gsub("[[:punct:]]", "-", spe))) # delete all special character from species name
        col_date <- stringr::str_which(names(dfy), "DATE")
        year <- dfy[1, col_date]

        # rasterisation
        # rp <- rasterize(y, ref_ras, "NOM_FRANCA")
        rp <- rasterize(y, ref_ras, 1)
        names(rp) <- spe

        print(y$NOM_FRANCA)

        # send the raster to s3
        path_s3 <- paste0("s3://bq-io/acer/melccfp/", group, "/melccfp_", spe, "_", year, ".tif")
        path <- paste0("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/rasterisation/melccfp_", spe, "_", year, ".tif")


        writeRaster(rp,
            path,
            filetype = "GTiff",
            gdal = "COMPRESS=DEFLATE",
            overwrite = TRUE
        )
        bash <- paste("AWS_ACCESS_KEY_ID=NJBPPQZX7PFUBP1LH8B0 AWS_SECRET_ACCESS_KEY=DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd cp -acl public-read", path, path_s3, sep = " ")
        system(bash)
        system(paste("rm", path))
    })
})

#### Creation des cartes de RS par groupe/sous-groupe(?) taxonomique
# retrieve the species names
system("AWS_ACCESS_KEY_ID=NJBPPQZX7PFUBP1LH8B0 AWS_SECRET_ACCESS_KEY=DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd ls s3://bq-io/acer/melccfp/reptiles/ > melccfp_reptiles_species.txt")
system("AWS_ACCESS_KEY_ID=NJBPPQZX7PFUBP1LH8B0 AWS_SECRET_ACCESS_KEY=DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd ls s3://bq-io/acer/melccfp/mammiferes/ > melccfp_mammiferes_species.txt")
system("AWS_ACCESS_KEY_ID=NJBPPQZX7PFUBP1LH8B0 AWS_SECRET_ACCESS_KEY=DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd ls s3://bq-io/acer/melccfp/poissons/ > melccfp_poissons_species.txt")
system("AWS_ACCESS_KEY_ID=NJBPPQZX7PFUBP1LH8B0 AWS_SECRET_ACCESS_KEY=DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd ls s3://bq-io/acer/melccfp/amphibiens/ > melccfp_amphibiens_species.txt")

spe_poiss <- read.table("melccfp_poissons_species.txt", h = F)[, 4]
spe_mamm <- read.table("melccfp_mammiferes_species.txt", h = F)[, 4]
spe_rept <- read.table("melccfp_reptiles_species.txt", h = F)[, 4]
spe_amph <- read.table("melccfp_amphibiens_species.txt", h = T)[, 4]

spe_list <- list(spe_poiss, spe_mamm, spe_rept, spe_amph)
group <- c("poissons", "mammiferes", "reptiles", "amphibiens")

for (l in 1:length(spe_list)) {
    spe_path <- c()
    x <- spe_list[[l]]
    group_name <- group[l]

    for (i in 1:length(x)) {
        path <- paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/melccfp/", group_name, "/", x[i])
        spe_path <- c(spe_path, path)
    }

    print(paste0("-----> Processing start for ", group_name))
    ras <- rast(spe_path)
    rs <- sum(ras, na.rm = T)

    # send the raster to s3
    path_s3 <- paste0("s3://bq-io/acer/ebv/rs_melccfp_", group_name, ".tif")
    path <- paste0("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/rasterisation/rs_melccfp_", group_name, ".tif")


    writeRaster(rs,
        path,
        filetype = "GTiff",
        gdal = "COMPRESS=DEFLATE",
        overwrite = TRUE
    )
    bash <- paste("AWS_ACCESS_KEY_ID=NJBPPQZX7PFUBP1LH8B0 AWS_SECRET_ACCESS_KEY=DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc S3_ENDPOINT_URL=https://object-arbutus.cloud.computecanada.ca s5cmd cp -acl public-read", path, path_s3, sep = " ")
    system(bash)
    system(paste("rm", path))
}

### test

rast <- rast("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/ebv/rs_melccfp_amphibiens.tif")
mapview::mapview(rast)
