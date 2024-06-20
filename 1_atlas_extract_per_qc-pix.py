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
