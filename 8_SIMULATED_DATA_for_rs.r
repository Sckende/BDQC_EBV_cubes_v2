# Obj : simulation of distribution map to dev script for the comparaison btw RS maps

library(mvtnorm) ## for dmvnorm()
library(viridis)
library(knitr)
library(animation)
ani.options(interval = 0.2, nmax = 200, autobrowse = FALSE) # animation aptions

#' @param coords a matrix-like structure of (x,y) coordinates
#' @param params list of parameters
simDens <- function(coords, params) {
    with(
        params,
        N * dmvnorm(x = coords, mean = mu, sigma = Sigma)
    )
}
#' @param time numeric time value
#' @param hyperpars list of hyper-parameters describing instantaneous config
paramFun <- function(time, hyperpars) {
    with(
        hyperpars,
        list(
            N = N, # nombre d'observations
            mu = c(
                x0 + ((time - xt0) / x1)^2, # moyenne pour x
                y0 + ((time - yt0) / y1)^2 # moyenne pour y
            ),
            Sigma = outer(c(sdx, sdy), c(sdx, sdy)) * matrix(c(1, rxy, rxy, 1), 2) # matrice de covariance
        )
    )
}
# -----> Exemple Guillaume
# parametres
dx <- dy <- 0.1
xmin <- ymin <- 0
xmax <- ymax <- 10
xvec <- seq(xmin, xmax, by = dx)
yvec <- seq(ymin, ymax, by = dy)
coords <- expand.grid(x = xvec, y = yvec)
tvec <- seq(0, 20, by = 0.1)
nx <- length(xvec)
ny <- length(yvec)

hpars <- list(
    x0 = 0,
    y0 = 0,
    x1 = 2,
    y1 = 1,
    xt0 = 5,
    yt0 = 2,
    sdx = 0.5, # taille blob
    sdy = 0.2, # taille blob
    rxy = 0.5,
    N = xmax * ymax
)

p0 <- paramFun(1, hpars) # parametres pour la simulation des valeurs dans la matrice
xx <- simDens(coords, p0) # simulation valeurs contenues dans la matrice
image(matrix(xx, nx, ny),
    col = viridis_pal()(20),
    useRaster = TRUE
)

dd <- array(
    dim = c(length(tvec), nx, ny),
    dimnames = list(t = tvec, x = xvec, y = yvec)
)
for (i in seq_along(tvec)) {
    dd[i, , ] <- matrix(simDens(coords, paramFun(tvec[i], hpars)), nx, ny)
}


saveHTML(
    {
        for (i in seq_along(tvec)) {
            image(dd[i, , ],
                col = viridis_pal()(20), useRaster = TRUE
            )
        }
    },
    movie.name = "moveblob.html"
)

# -----> Application aire de distri
# raster en exemple
library(terra)
espece <- "setophaga_virens"
# --> choix methode
# methode <- "Maxent"
# methode <- "brt"
methode <- "randomForest"
# chargement des rasters
ras_ex <- rast(paste0("https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/TdB_bench_maps/maps/CROPPED_QC_", espece, "_", methode, "_Predictors_Bias_Spatial.tif"))
# xmin = -832462.4
# xmax = 785810
# ymin = 120547.3
# ymax = 2088446

# parametres
# dx1 <- dy1 <- res(ras_ex)[1]
# xmin1 <- ext(ras_ex)[1] # for coord
# ymin1 <- ext(ras_ex)[3] # for coord
# xmax1 <- ext(ras_ex)[2] # for coord
# ymax1 <- ext(ras_ex)[4] # for coord
dx1 <- dy1 <- 1
xmin1 <- 1 # for coord
ymin1 <- 1 # for coord
xmax1 <- 162 # for coord
ymax1 <- 197 # for coord
xvec1 <- seq(xmin1, xmax1, by = dx1) # coord - longitudes
yvec1 <- seq(ymin1, ymax1, by = dy1) # coord - latitudes
coords1 <- expand.grid(x = xvec1, y = yvec1)
tvec1 <- seq(1, 25, by = 1) # time - third dimension of the matrix
nx1 <- length(xvec1)
ny1 <- length(yvec1)

hpars1 <- list(
    x0 = 0,
    y0 = 0,
    x1 = 2,
    y1 = 1,
    xt0 = 5,
    yt0 = 2,
    sdx = 100, # taille blob
    sdy = 100, # taille blob
    rxy = 0.5,
    N = xmax1 * ymax1
)

p01 <- paramFun(1, hpars1) # parametres pour la simulation des valeurs dans la matrice
xx1 <- simDens(coords1, p01) # simulation valeurs contenues dans la matrice
image(rast(matrix(xx1, nx1, ny1)),
    col = viridis_pal()(20),
    useRaster = TRUE
)

dd1 <- array(
    dim = c(length(tvec1), nx1, ny1),
    dimnames = list(t = tvec1, x = xvec1, y = yvec1)
)
for (i in seq_along(tvec1)) {
    dd1[i, , ] <- matrix(simDens(coords1, paramFun(tvec1[i], hpars1)), nx1, ny1)
}

x11()
par(mfrow = c(5, 5))
for (i in seq_along(tvec1)) {
    plot(rast(dd1[i, , ]))
}

# saveHTML(
#     {
#         for (i in seq_along(tvec1)) {
#             image(rast(dd1[i, , ]),
#                 col = viridis_pal()(20), useRaster = TRUE
#             )
#         }
#     },
#     movie.name = "moveblob.html"
# )


# ----- > autre technique ?
library(stpp)
sim.stpp(class = "poisson", s.region = coords, npoints = 1000, nsim = 1)


#### generation de cartes avec maxent ####

# test avec Antigone canadensis

# library(ratlas)
library(dplyr)
library(duckdb)
library(duckdbfs)
library(dbplyr)
library(sf)
library(mapview)
library(terra)

source("http://atlas.biodiversite-quebec.ca/bq-atlas-parquet.R")
atlas_dates
atlas_rem <- atlas_remote(tail(atlas_dates$dates, n = 1))

spe <- c("Aix sponsa", "Branta canadensis", "Anas rubripes", "Anas acuta", "Anas platyrhynchos", "Spatula discors")

spe_obs <- atlas_rem |>
    filter(valid_scientific_name %in% spe) |>
    select(valid_scientific_name, latitude, longitude, dataset_name, year_obs, month_obs, day_obs) |>
    collect()

dim(spe_obs)
colnames(spe_obs)

spe_obs |>
    group_by(valid_scientific_name, year_obs) |>
    arrange(year_obs) |>
    summarize(cnt = n()) |>
    print(n = 150)

spe_obs2 <- spe_obs[spe_obs$year_obs %in% 1990:2020 & spe_obs$month_obs %in% 4:7, ]


# Conversion to sf object
spe_obs_sf <- st_as_sf(spe_obs2,
    coords = c("longitude", "latitude"),
    crs = 4326
)
mapview(split(spe_obs_sf, spe_obs_sf$valid_scientific_name))
