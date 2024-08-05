# Objectifs: produire des fonctions pour extraire l'information de chaque axe du cube
# taxonomique
# spatial
# twemporel
# Exercice effectue sur les cartes de Vincent Bellavance

library(terra)

# test avec une seule map
path <- "https://object-arbutus.cloud.computecanada.ca/bq-io/acer/oiseaux-nicheurs-qc/melospiza_lincolnii_range_1999.tif"

mp <- rast(path)
mp

null_ras <- mp
values(null_ras) <- 0

# recupere les ids de pixels avec une valeur = 1
ids <- cells(mp, 1)[[1]]
ids

df <- data.frame(
    sp = "melospiza_lincolnii",
    year = 1999,
    id_pix = ids,
    spe_rich = 1
)

null_ras[as.integer(df$id_pix)] <- df$spe_rich
plot(null_ras)


#### Exploration du package ebvcube ####
# https://cran.r-project.org/web/packages/ebvcube/readme/README.html
library("ebvcube")

# set the path to the file
file <- system.file(file.path("extdata", "martins_comcom_subset.nc"), package = "ebvcube")

# read the properties of the file
prop.file <- ebv_properties(file, verbose = FALSE)

# take a look at the general properties of the data set - there are more properties to discover!
prop.file@general[1:4]

# datacube associated with the netcdf file
datacubes <- ebv_datacubepaths(file, verbose = FALSE)
datacubes
# In the next step we will get the properties of one specific datacube - fyi: the result also holds the general file properties from above.
prop.dc <- ebv_properties(file, datacubes[1, 1], verbose = FALSE)
prop.dc@metric

# plot the global map
dc <- datacubes[2, 1]
ebv_map(file, dc,
    entity = 1, timestep = 12, classes = 9,
    verbose = FALSE, col_rev = TRUE
)

# It’s nice to see the global distribution, but how is the change of that datacube (non forest birds) over time? Let’s take a look at the average. The function returns the values, catch them!

# get the averages and plot
averages <- ebv_trend(file, dc, entity = 1, verbose = FALSE)
averages

# It would be cool to have that for other indicators as well? Check out the different options for ‘method’.
# Before you actually load the data it may be nice to get an impression of the value range and other basic measurements.

# info for whole dataset
measurements <- ebv_analyse(file, dc, entity = 1, verbose = FALSE)
# see the included measurements
names(measurements)

# check out the mean and the number of pixels
measurements$mean
#> [1] 0.6140731
measurements$n
#> [1] 7650

# info for a subset defined by a bounding box
# you can also define the subset by a Shapefile - check it out!
bb <- c(-26, 64, 30, 38)
measurements.bb <- ebv_analyse(file, dc, entity = 1, subset = bb, verbose = FALSE)
# check out the mean of the subset
measurements.bb$mean
#> [1] 0.3241093
measurements.bb$n
#> [1] 720

# To access the first three timesteps of the data you can use the following:

# load whole data as array for two timesteps
data <- ebv_read(file, dc, entity = 1, timestep = 1:3, type = "a")
dim(data)
#> [1] 85 90  3

# You can also get a spatial subset of the data by providing a Shapefile.

# load subset from shapefile (Cameroon)
shp <- system.file(file.path("extdata", "cameroon.shp"), package = "ebvcube")
data.shp <- ebv_read_shp(file, dc, entity = 1, shp = shp, timestep = c(1, 2, 3), verbose = FALSE)
dim(data.shp)

# very quick plot of the resulting raster plus the shapefile
borders <- terra::vect(shp)
ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = data.shp[[1]]) +
    tidyterra::geom_spatvector(data = borders, fill = NA) +
    ggplot2::scale_fill_fermenter(na.value = NA, palette = "YlGn", direction = 1) +
    ggplot2::theme_classic()
