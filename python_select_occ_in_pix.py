import pandas as pd
import geopandas as geopd
import numpy as np
import matplotlib.pyplot as plt

# importation du parquet local via pandas
atlas = pd.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_2024-05-29.parquet", engine = "pyarrow")
atlas.shape

atlas = atlas[atlas["within_quebec"] == "t"]
type(atlas)
atlas.dtypes
atlas.shape
print(atlas.head())

# Conversion du dataframe en geodataframe pour manipulation avec geopandas
geoatlas = geopd.GeoDataFrame(atlas, geometry = geopd.points_from_xy(atlas.longitude, atlas.latitude), crs = "EPSG:4326")
type(geoatlas)
geoatlas.shape

# verification with a sample
samp = geoatlas[0:10001]
samp.tail()
samp.shape
samp.plot()
plt.show()
samp.explore()


# importation du geopackage avec le Qc pixelisé
pix = geopd.read_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale.gpkg")

# test with teh first pixel of Qc
type(pix)
pix.shape
pix.head()
pix1 = pix.iloc[[0]]
type(pix1)
pix1
# conversion of the crs pix1
pix1 = pix1.to_crs(4326)
test_within = geopd.sjoin(geoatlas, pix1, op = "within")
test_within