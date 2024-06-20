# import pandas as pd
# import geopandas as geopd
# import numpy as np
# import matplotlib.pyplot as plt
# import geoparquet as gpq
import duckdb as d

# from shapely.ops import unary_union

# con = d.connect()
# con.install_extension("spatial")
# con.load_extension("spatial")

# avec duckdb

# Load spatial extension
d.sql("INSTALL spatial; LOAD spatial")
# Load HTTPs extension for request on s3 buckets
d.sql("INSTALL https; LOAD https")

# Lecture Atlas (parquet)
atlas = d.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_2024-05-29.parquet")

#https://object-arbutus.cloud.computecanada.ca/bq-io/acer/TdeB_benchmark_SDM/QUEBEC_Unique_poly.gpk

d.sql("CREATE SECRET secret1 (TYPE S3, KEY_ID 'NJBPPQZX7PFUBP1LH8B0', SECRET 'DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc', ENDPOINT 's3.object-arbutus.cloud.computecanada.ca')")
d.sql("DROP SECRET secret1")

d.sql("CREATE TABLE atlas AS SELECT * FROM 's3://bq-io/atlas/parquet/atlas_2024-05-29.parquet'")


# d.sql("CREATE TABLE atlas_sf AS SELECT *, ST_Point(CAST(longitude as float), CAST(latitude as float)) AS geometry, FROM atlas")

# # Lecture qc-pix (geopackage)
# d.sql("CREATE TABLE qc_pix AS SELECT * FROM ST_Read('/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale_latlon.gpkg')")

# # Merge Atlas & qc-pix with a loop
# # --> exploration
# d.sql("SELECT min(year_obs) FROM atlas_sf") # 1534
# d.sql("SELECT max(year_obs) FROM atlas_sf") #2023

# d.sql("SELECT year_obs, COUNT(year_obs) FROM atlas_sf WHERE CAST(year_obs AS int) >= 1900 GROUP BY year_obs ORDER BY year_obs DESC")

# # --> loop creation
# for year in range(1900, 2024):
#     print("---------------------------------------------->" + str(year))
#     qu = "CREATE TABLE atlas_pix AS SELECT * FROM atlas_sf, qc_pix WHERE atlas_sf.year_obs == " + str(year) + " AND ST_Intersects(atlas_sf.geometry, qc_pix.geom)"
#     print(qu)
#     d.sql(qu)
#     # Write the new db
#     qu2 = "COPY (SELECT * FROM atlas_pix) TO 'atlas_pix_parquet/atlas_pix_" + str(year) + ".parquet' (FORMAT 'parquet')"
#     d.sql(qu2)
#     d.sql("DROP TABLE atlas_pix")






























# gestion des tables
d.sql("DROP TABLE atlas_sf")


d.sql("SELECT DISTINCT year_obs FROM atlas_pix")

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