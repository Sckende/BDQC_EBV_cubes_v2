import duckdb as d

# Load spatial extension
d.sql("INSTALL spatial; LOAD spatial")
# Load HTTPs extension for request on s3 buckets
d.sql("INSTALL https; LOAD https")

# Lecture Atlas (parquet) local 
atlas = d.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_2024-05-29.parquet")

# Lecture Atlas en remote
d.sql("CREATE SECRET secret1 (TYPE S3, KEY_ID 'XXXX', SECRET 'xxxxx', URL_STYLE path, ENDPOINT 'object-arbutus.cloud.computecanada.ca')")
d.sql("CREATE TABLE atlas AS SELECT * FROM 's3://bq-io/atlas/parquet/atlas_2024-05-29.parquet'")
# d.sql("DROP SECRET secret1")

# Conversion en objet spatial
d.sql("CREATE TABLE atlas_sf AS SELECT *, ST_Point(CAST(longitude as float), CAST(latitude as float)) AS geometry, FROM atlas")

# Lecture qc-pix (geopackage)
d.sql("CREATE TABLE qc_pix AS SELECT * FROM ST_Read('/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale_latlon.gpkg')")

# Merge Atlas & qc-pix with a loop
# --> exploration
d.sql("SELECT min(year_obs) FROM atlas_sf") # 1534
d.sql("SELECT max(year_obs) FROM atlas_sf") #2023

d.sql("SELECT year_obs, COUNT(year_obs) FROM atlas_sf WHERE CAST(year_obs AS int) >= 1900 GROUP BY year_obs ORDER BY year_obs DESC")

# Trying to get all column names
d.sql("DESCRIBE atlas_sf LIMIT 28")
d.sql("SHOW atlas_sf")
d.sql("SELECT COLUMN_NAME from information_schema.columns WHERE table_name = 'atlas_sf'")
d.sql("SELECT * FROM atlas_sf").show(max_width=280)

d.sql("SHOW qc_pix")

# --> loop creation
for year in range(1900, 2024):
    print("---------------------------------------------->" + str(year))
    qu = "CREATE TABLE atlas_pix AS SELECT atlas_sf.valid_scientific_name, atlas_sf.longitude, atlas_sf.latitude, atlas_sf.year_obs, atlas_sf.month_obs, atlas_sf.day_obs, atlas_sf.group_en, atlas_sf.kingdom, atlas_sf.phylum, atlas_sf.class, atlas_sf.order, atlas_sf.family, atlas_sf.genus, qc_pix.ID AS pix_id, qc_pix.geom AS geom_sql, to_binary(qc_pix.geom) AS geometry FROM atlas_sf, qc_pix WHERE atlas_sf.year_obs == " + str(year) + " AND ST_Intersects(atlas_sf.geometry, qc_pix.geom)"
    print(qu)
    d.sql(qu)
    # d.sql("SHOW atlas_pix")
    # d.sql("SELECT * FROM atlas_pix")


    # Write the new db
    qu2 = "COPY (SELECT * FROM atlas_pix) TO 'atlas_pix_parquet/atlas_pix_" + str(year) + ".parquet' (FORMAT 'parquet')"
    d.sql(qu2)
    d.sql("DROP TABLE atlas_pix")








### WARNING !!! ###
# voir avec Guillaume pour la création d'un géoparquet valide, c'est à dire créer le geo metadata
# voir également pourquoi impossible de convertir le blob créé à partir de geomn dans SQL dans python 





# ------------------------------------------------------------- #
# Zone de verification
year=2015
qu = "CREATE TABLE atlas_pix AS SELECT * FROM atlas_sf, qc_pix WHERE atlas_sf.year_obs == " + str(year) + " AND ST_Intersects(atlas_sf.geometry, qc_pix.geom)"
print(qu)
d.sql(qu)

d.sql("SELECT * FROM atlas_pix")



# Zone de test pour writer un file en geoparquet
# Etapes de Guillaume: requete sql sur db, write in csv, converti le csv en parquet avec GDAL
# cf script geospatial-docker/atlas2pmtiles.sh 
d.sql("COPY (SELECT * FROM atlas_pix) TO 'atlas2015_test.csv' CSV HEADER DELIMITER ','")
#GDAL
ogr2ogr -f Parquet -s_srs EPSG:4326 -t_srs EPSG:4326 atlas2015_test.parquet atlas2015_test.csv

#--> voir avec Guillaume
# pb est qu'il manque REQUIRED METADATA pour valider le parquet en geoparquet : A GeoParquet file MUST include a geo key in the Parquet metadata (see FileMetaData::key_value_metadata). The value of this key MUST be a JSON-encoded UTF-8 string representing the file and column metadata that validates against the GeoParquet metadata schema. The file and column metadata fields are described below (cf https://geoparquet.org/releases/v1.0.0/ & https://github.com/geopandas/geopandas/discussions/3158) 
import geopandas as geopd
import pandas as pd

test2 = geopd.read_parquet("atlas2015_test.parquet")
at2015 = pd.read_parquet("atlas2015_test.parquet")
at2015.columns.values
pix2015 = at2015[['ID', 'geom']]
from shapely import wkt
pix2015["geom"] = geopd.GeoSeries.from_wkt(pix2015["geom"])
pix2015_sf = geopd.GeoDataFrame(pix2015, geometry = 'geom', crs = "EPSG:4326")
type(pix2015_sf)

pix2015_sf.iloc[[1]].plot()
pix2015_sf.plot()
plt.show()
type(pix2015_sf.iloc[[1]])

# ----- #
occ2015 = at2015.iloc[:,0:27]
type(occ2015)
occ2015.columns.values
occ2015["geometry"] = geopd.GeoSeries.from_wkt(occ2015["geometry"])
occ2015_sf = geopd.GeoDataFrame(occ2015, geometry = 'geometry', crs = "EPSG:4326")
occ2015_sf.shape
occ2015_sf.plot()
plt.show()
# ID conversion to numeric
occ2015_sf.ID = pd.to_numeric(occ2015_sf.ID)

occ2015_sf.ID.unique()
occ2015_sf.ID.value_counts().sort_index()
occ2015_sf.geometry.value_counts().sort_index()

max(occ2015_sf.ID.value_counts())
occ2015_sf[occ2015_sf['ID'] == 26].shape
pix26=occ2015_sf[occ2015_sf['ID'] == 26]
pix26.latitude.describe()
pix26.longitude.describe()

pix26[['longitude', 'latitude', 'geometry']]

import matplotlib.pyplot as plt
pix26.plot()
plt.show()

occ2015_sf[occ2015_sf['ID'] == '44840'].plot()
pix2015_sf.iloc[[1]].plot()

plt.show()


#### Test from the generated parquet
import geopandas as geopd
import pandas as pd
import shapely as shp
import matplotlib.pyplot as plt

pq1984=pd.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/atlas_pix_parquet/atlas_pix_1984.parquet")
type(pq1984.geom)
pq1984.columns.values
pq1984.dtypes
pq1984.geom.describe()
pq1984.geometry.describe()
pq1984.geometry_text.describe()
pq1984[['geometry', 'geom']]

test=pq1984[['geometry_text']]
test_sf=shp.from_wkb(test)


test_sf.plot()
plt.show()
test.describe()

# Erreur incomprehensible
# voir ressources https://shapely.readthedocs.io/en/stable/reference/shapely.from_wkb.html
t=pq1984.geometry[0]
shp.from_wkb(t)
pq1984["geometry_bin"] = geopd.GeoSeries.from_wkb(pq1984["geometry_bin"])
pq1984["geometry"] = geopd.GeoSeries.from_wkb(pq1984["geometry"])

pq1984.geometry_bin.describe()

pq1984_sf = geopd.GeoDataFrame(pq1984, geometry = 'geom', crs = "EPSG:4326")
pq1984_sf.geom.describe()
pq1984_sf.plot()
plt.show()



test2 = geopd.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/atlas_pix_parquet/atlas_pix_1984.parquet")
