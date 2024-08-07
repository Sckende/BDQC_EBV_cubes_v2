import pandas as pd
import geopandas as geopd
import numpy as np
import matplotlib.pyplot as plt
import geoparquet as gpq
import duckdb as d

#### Explo duckDB ####
# gestion des tables
d.sql("DROP TABLE atlas_sf")

d.sql("SELECT DISTINCT year_obs FROM atlas")
d.sql("SELECT * FROM atlas")


# test
test = d.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/atlas_pix_parquet/atlas_pix_2011.parquet")

d.sql("SELECT geom FROM test")

data2000 = d.sql("SELECT * FROM dataB WHERE year_obs == '2000'")
d.sql("SELECT CAST(latitude as float) FROM data2000")

d.sql("CREATE TABLE data2000_sf AS SELECT *, ST_Point(CAST(longitude as float), CAST(latitude as float)) AS geometry, FROM data2000")
d.sql("SELECT * FROM data2000_sf")

d.sql("CREATE TABLE data2018_sf AS SELECT *, ST_Point(CAST(longitude as float), CAST(latitude as float)) AS geometry, FROM dataB WHERE year_obs == '2018'")
d.sql("SELECT * FROM data2018_sf")

# Lecture qc-pix
d.sql("CREATE TABLE multi_pix AS SELECT * FROM ST_Read('/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale_latlon.gpkg')")
d.sql("SELECT COUNT(ID) FROM multi_pix")

# st_intersect between pix and occs 2000
d.sql("SELECT data2000_sf.valid_scientific_name FROM data2000_sf, pix2 WHERE ST_Intersects(data2000_sf.geometry, pix2.geom)")

# st_intersect between pix and occs 2018
d.sql("SELECT data2018_sf.valid_scientific_name FROM data2018_sf, pix3 WHERE ST_Intersects(data2018_sf.geometry, pix3.geom)")

# st_intersect between multipix and occs 2000
d.sql("CREATE TABLE union2000 AS SELECT * FROM data2000_sf, multi_pix WHERE ST_Intersects(data2000_sf.geometry, multi_pix.geom)")
d.sql("SELECT COUNT(DISTINCT family) FROM union2000 WHERE ID == 965405")

# Extra de gestion
d.sql("SHOW TABLES")
d.sql("DROP TABLE multi_pix")

d.sql("DESCRIBE data2000_sf")
d.sql("SELECT COLUMN_NAME from information_schema.columns WHERE table_name = 'union2000'")

d.sql("INSTALL spatial; LOAD spatial")
d.sql("CREATE TABLE test_sf AS SELECT *, ST_asWKB(ST_Point(longitude, latitude)) as geometry, FROM data2000;")

#### Exploration pandas, geopandas, matplotlib ####

# --> Union Atlas-QcPix
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
geoatlas.dtypes
pd.Series.value_counts(pd.isna(geoatlas.within_quebec))

# Conversion & storage in GEOPARQUET
geoatlas.to_geoparquet("atlas.geoparquet")
# pre-treatmentr for excluding points outside of Qc
qc = geopd.read_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/QUEBEC_CR_NIV_01.gpkg").to_crs(4326)
qc
qc_fus = geopd.GeoDataFrame(geopd.GeoSeries(unary_union(qc.geometry)))
qc_fus = qc_fus.rename(columns={0:"geometry"}).set_geometry("geometry").set_crs(4326)
qc_fus.plot()
plt.show()

base = qc.plot()
test.plot(ax = base, color = "grey")
plt.show()
# verification with a sample
samp = geoatlas[0:1000001]
samp.tail()
samp.shape
samp.plot()
plt.show()
samp.explore()

samp.to_geoparquet("samp.geoparquet")



# importation du geopackage avec le Qc pixelisé
pix = geopd.read_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale.gpkg")

# test with teh first pixel of Qc
type(pix)
pix.shape
pix.head()

# conversion of the crs pix
pix = pix.to_crs(4326)

# test extraction point
# pix1 = pix.iloc[[102]]
pix1 = pix[0:51]
type(pix1)
pix1

# visualization
base = samp.plot()
pix1.plot(ax = base, edgecolor = "red")
plt.show()

# extraction test with small sample
test_within = geopd.sjoin(samp, pix1)
test_within
test_within.dtypes
test_within.index_right #les deux dernieres colonnes correspondent à l'ID du point et à l'ID du polygone1
test_within.shape
print(test_within)

# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
#     print(test_within)

# extraction test with small smple of poly-pixel and atlas in 2000
# car crash lors du test sur atlas au complet
geoatlas.dtypes
gatlas2000 = geoatlas[pd.to_numeric(geoatlas.year_obs) == 2000]
test = geopd.sjoin(gatlas2000, qc_fus) # treatment - faster with gdal ?
gatlas2000
gatlas2000.columns.values
gatlas2000.shape
gatlas2000.year_obs.value_counts()
test_within = geopd.sjoin(gatlas2000, pix1, how = "left") # CRASH sur atlas au complet !
test_within.shape
# with pd.option_context(10, None, 'display.max_columns', None):  # more options can be specified also
    # print(test_within)

pd.set_option("display.max_columns", 28)
print(test_within)
test_within.ID.value_counts()
nann = test_within[pd.isna(test_within["ID"])]
nann.shape

nona = test_within.dropna()
nona.shape
nona.ID.value_counts()

base = gatlas2000.plot(color = "orange")
nona.plot(ax = base, color = "red")
plt.show()


### Sandbox pour confirmer la méthode --> qui est confirmée !
import geodatasets

chicago = geopd.read_file(

    geodatasets.get_path("geoda.chicago_health")

)

groceries = geopd.read_file(

    geodatasets.get_path("geoda.groceries")

).to_crs(chicago.crs)

groceries_w_communities = geopd.sjoin(groceries, chicago)
groceries_w_communities
groceries_w_communities.columns.values
groceries_w_communities.ComAreaID.unique()
groceries_w_communities.ComAreaID.describe()
groceries_w_communities.ComAreaID.value_counts()



groceries_w_communities.head()
chicago.shape  
groceries.shape
groceries_w_communities.shape

base = chicago.plot()
groceries.plot(ax = base, color = "red")
chicago[chicago["ComAreaID"].isin(groceries_w_communities.ComAreaID.unique())].plot(ax = base, color = "green")
plt.show()

print(groceries_w_communities.dtypes)

## --> Lecture & manipulation des fichiers parquets QcPix
atlas2015 = pd.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/atlas_pix_parquet/atlas_pix_2015.parquet", engine = "pyarrow")
atlas2015.columns.values
type(atlas2015)
atlas2015.geometry

at2015 = geopd.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/atlas_pix_parquet/atlas_pix_2015.parquet")

#### Test from the generated parquet ####
import geopandas as geopd
import pandas as pd
import shapely as shp
import matplotlib.pyplot as plt

pq1984=pd.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_pix_parquet/atlas_pix_2020.parquet")
pq1984.columns.values
type(pq1984.geometry)
pq1984.dtypes
pq1984.geometry.describe()
pq1984.pix_id.describe()

# spatial object conversion

# ---> keep the unique pixels
uniq_pix = pq1984.drop_duplicates("geometry")
uniq_pix.columns.values
uniq_pix.geometry
type(uniq_pix)
uniq_pix.dtypes
uniq_pix.shape
uniq_pix["geometry"] = geopd.GeoSeries.from_wkt(uniq_pix["geometry"]) # Voir avec Guillaume pourquoi cette ligne est indispensable?
uniq_pix_sf = geopd.GeoDataFrame(uniq_pix, geometry = "geometry", crs = "EPSG:4326")
type(uniq_pix_sf)
uniq_pix_sf.shape
pix1984=uniq_pix_sf[["pix_id", "geometry"]]
type(pix1984)
pix1984.dtypes

# ---> visualisation
pix1984.plot()
plt.show()

# ---> Compute the species richness by pixel
# subtest = pq1984.iloc[1:100]
# subtest.geometry.describe()
# subtest.columns.values

l=pq1984.groupby("pix_id")["valid_scientific_name"].unique()
info=pq1984.groupby("pix_id")["valid_scientific_name"].nunique()
type(info)
# conversion from series to dataframe
info_df = pd.DataFrame({'pix_id':info.index, 'spe_rich':info.values, 'spe_list':l.values})
type(info_df)
# retrieve the geom pixel
# pix_geom=subtest[["pix_id", "geometry"]]
# type(pix_geom)
# res=pd.merge(left=info_df, right=pix_geom, on="pix_id", how="left")
# res.columns.values

# pix1984.dtypes
# info_df.dtypes
# info_df.join(pix1984, on="pix_id", how="left")

info_df2=info_df.merge(pix1984, how='left', on='pix_id') #dataframe
# conversion to geodataframe
type(info_df2)
info_sf=geopd.GeoDataFrame(info_df2, geometry = "geometry", crs = "EPSG:4326")
# ---> visualisation - choropleth map (maps where the color of each shape is based on the value of an associated variable)
#https://geopandas.org/en/stable/docs/user_guide/mapping.html

info_sf.plot(column="spe_rich", legend=True, cmap="OrRd", scheme="quantiles")
plt.show()


#### Travail sur les rasters ####

import geopandas as gpd
import matplotlib.pyplot as plt
# from matplotlib import pyplot
import numpy as np
import rasterio as rio
from rasterio import plot
from rasterio.plot import show
from rasterio.mask import mask

# get polygons
shapef=gpd.read_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_regions/CERQ_SHP/CR_NIV_01_S.shp")
shapef.plot(column="FID01")
plt.show()

# get raster
raster1=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/vireo_flavifrons_range_2017.tif")
plt.imshow(raster1.read(1))
plt.show()
raster1.meta


# get the crs
shapef.crs #321988
raster1.crs #32198

# convert crs polygons to crs raster
sf2=shapef.to_crs(32198)
sf2.plot(facecolor="none", edgecolor="grey")
plt.show()
 
# plot raster & polygons side by side
fig, (axr, axg) = plt.subplots(1, 2)
show((raster1, 1), ax=axr, cmap='Reds', title="distribution")
sf2.plot(ax=axg, facecolor="none", edgecolor="grey")
plt.show()

# plot polygons over raster
fig, (ax1, ax2) = plt.subplots(1,2)
show((raster1, 1), ax=ax1, cmap='Reds', title="distribution")
sf2.plot(ax=ax1, facecolor="none", edgecolor="grey", linewidth = 0.25)
sf2.iloc[[2]].plot(ax=ax2)
plt.show()

# ---> look https://www.youtube.com/watch?v=Tqph7_qMujk & https://www.youtube.com/watch?v=LRbOgRKWVag for subplot informationsfig, (ax1, ax2) = plt.subplots(1,2)
show((sp1, 1), ax=ax1)
show((sp2, 1), ax=ax2)
plt.show()

#### clip raster from a polygon ####
# -------------------------------- # 
import pycrs
# Read the polygon
poly=sf2.iloc[[3]]
type(poly)

# read the raster
raster1

# Check if the crs is the same for the two files
poly.crs
raster1.crs

# crop
rastmask,_=rio.mask.mask(raster1, poly.geometry, crop=True)

# Visualisation
fig = plt.figure(figsize=[12,8])
# Plot the raster data using matplotlib
ax = fig.add_axes([0, 0, 1, 1])
raster_image=ax.imshow(rastmask[0,:,:])
plt.show()

# delete nan values
type(rastmask)
rast2=rastmask[~np.isnan(rastmask)]
rast2.sum()
rast2.max()

#### Stack several rasters and add description in metadata (species names) ####
#1-Modify the metadata to include the species names
raster1=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/vireo_flavifrons_range_2017.tif")
metadata=raster1.meta.copy()
type(metadata)
metadata['name'] = "Species name"

raster1.meta.update() = metadata
raster1.tags()


#### Stack several rasters ####
# -------------------------- # 
import rasterio as rio
import matplotlib.pyplot as plt
from rasterio.plot import show


sp1=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/vireo_flavifrons_range_2017.tif")
plt.imshow(sp1.read(1))
plt.show()
sp2=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/vireo_gilvus_range_2017.tif")
sp3=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/vireo_olivaceus_range_2017.tif")
sp4=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/vireo_philadelphicus_range_2017.tif")
sp5=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/vireo_solitarius_range_2017.tif")

# visual
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3)
show((sp1, 1), ax=ax1)
show((sp2, 1), ax=ax2)
show((sp3, 1), ax=ax3)
show((sp4, 1), ax=ax4)
show((sp5, 1), ax=ax5)
plt.show()

# Preparation pour le stack
out_img="stack_sp.tif"

out_meta=sp5.meta.copy()
out_meta.update({"count": 5})

# stackage
file_list=[sp1, sp2, sp3, sp4, sp5]
with rio.open(out_img, 'w', **out_meta) as dest:
    for band_nr, src in enumerate(file_list, start=1):
        dest.write(src, band_nr)


stack=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/stack_sp.tif")
type(stack)


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3)
show((stack, 1), ax=ax1)
show((stack, 2), ax=ax2)
show((stack, 3), ax=ax3)
show((stack, 4), ax=ax4)
show((stack, 5), ax=ax5)
plt.show()

show(stack.read(1))
plt.show()

### autre methode
from osgeo import gdal #impossible installer osgeo