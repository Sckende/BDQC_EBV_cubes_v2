import duckdb as d

# Load spatial extension
d.sql("INSTALL spatial; LOAD spatial")
# Load HTTPs extension for request on s3 buckets
d.sql("INSTALL https; LOAD https")

# Lecture Atlas (parquet) local 
atlas = d.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_2024-05-29.parquet")

# Lecture Atlas en remote
d.sql("CREATE SECRET secret1 (TYPE S3, KEY_ID 'NJBPPQZX7PFUBP1LH8B0', SECRET 'DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc', URL_STYLE path, ENDPOINT 'object-arbutus.cloud.computecanada.ca')")
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

# --> loop creation
for year in range(1900, 2024):
    print("---------------------------------------------->" + str(year))
    qu = "CREATE TABLE atlas_pix AS SELECT * FROM atlas_sf, qc_pix WHERE atlas_sf.year_obs == " + str(year) + " AND ST_Intersects(atlas_sf.geometry, qc_pix.geom)"
    print(qu)
    d.sql(qu)

    # Write the new db
    qu2 = "COPY (SELECT * FROM atlas_pix) TO 'atlas_pix_parquet/atlas_pix_" + str(year) + ".parquet' (FORMAT 'parquet')"
    d.sql(qu2)
    d.sql("DROP TABLE atlas_pix")
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
test2 = geopd.read_parquet("atlas2015_test.parquet")
at2015 = pd.read_parquet("atlas2015_test.parquet")
at2015.columns.values
pix2015 = at2015[['ID', 'geom']]
from shapely import wkt
pix2015["geom"] = geopd.GeoSeries.from_wkt(pix2015["geom"])
pix2015_sf = geopd.GeoDataFrame(pix2015, geometry = 'geom', crs = "EPSG:4326")

pix2015_sf.iloc[[1]].plot()
plt.show()
type(pix2015_sf.iloc[[1]])

# ----- #
occ2015 = at2015.iloc[:,0:27]
type(occ2015)
occ2015.columns.values
occ2015["geometry"] = geopd.GeoSeries.from_wkt(occ2015["geometry"])
occ2015_sf = geopd.GeoDataFrame(occ2015, geometry = 'geometry', crs = "EPSG:4326")
occ2015_sf.ID = pd.to_numeric(occ2015_sf.ID)
occ2015_sf.ID.unique()
occ2015_sf.ID.value_counts().sort_index()
max(occ2015_sf.ID.value_counts())
occ2015_sf[occ2015_sf['ID'] == '269495'].shape

occ2015_sf[occ2015_sf['ID'] == '44840'].plot()
pix2015_sf.iloc[[1]].plot()

plt.show()
