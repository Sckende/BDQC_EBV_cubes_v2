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
