import duckdb as d
import os



# Load spatial extension
d.sql("INSTALL spatial; LOAD spatial")
# Load HTTPs extension for request on s3 buckets
d.sql("INSTALL https; LOAD https")

# Lecture Atlas (parquet) local 
# atlas = d.read_parquet("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/atlas_2024-07-08.parquet")

# Lecture Atlas en remote
d.sql("CREATE SECRET secret1 (TYPE S3, KEY_ID 'NJBPPQZX7PFUBP1LH8B0', SECRET 'DVQZTIQYUBxqs0nwtfA4n1meL8Fv9w977pSp8Gjc', URL_STYLE path, ENDPOINT 'object-arbutus.cloud.computecanada.ca')")
# d.sql("CREATE TABLE atlas AS SELECT * FROM 's3://bq-io/atlas/parquet/atlas_2024-05-29.parquet'") # parquet file
d.sql("CREATE TABLE atlas_sf AS SELECT * FROM 's3://bq-io/atlas/parquet/atlas_2024-07-16.parquet'") # GEOparquet file
d.sql("SELECT DISTINCT year_obs FROM atlas_sf")
d.sql("SELECT group_en, COUNT(group_en) FROM atlas_sf GROUP BY group_en")
d.sql("SELECT group_fr, COUNT(group_fr) FROM atlas_sf GROUP BY group_fr")

# d.sql("DROP SECRET secret1")

# Conversion en objet spatial
d.sql("CREATE TABLE atlas_sf AS SELECT *, ST_Point(CAST(longitude as float), CAST(latitude as float)) AS geometry, FROM atlas")

# Lecture qc-pix (geopackage)
d.sql("CREATE TABLE qc_pix AS SELECT * FROM ST_Read('/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/qc_polygons/qc_grid_1x1km_finale_latlon.gpkg')")

# Merge Atlas & qc-pix with a loop
# --> exploration
d.sql("SELECT min(year_obs) FROM atlas_sf") # 1534
d.sql("SELECT max(year_obs) FROM atlas_sf") # 2023

d.sql("SELECT year_obs, COUNT(year_obs) FROM atlas_sf WHERE CAST(year_obs AS int) >= 1900 GROUP BY year_obs ORDER BY year_obs DESC")

# Trying to get all column names
d.sql("DESCRIBE atlas_sf LIMIT 28")
d.sql("SHOW atlas_sf")
d.sql("SELECT COLUMN_NAME from information_schema.columns WHERE table_name = 'atlas_sf'")
d.sql("SELECT * FROM atlas_sf").show(max_width=280)

d.sql("SHOW qc_pix")
d.sql("SHOW atlas_sf")

# --> loop creation
for year in range(1900, 1901):
    print("---------------------------------------------->" + str(year))
    qu = "CREATE TABLE atlas_pix AS SELECT atlas_sf.valid_scientific_name, atlas_sf.longitude, atlas_sf.latitude, atlas_sf.year_obs, atlas_sf.month_obs, atlas_sf.day_obs, atlas_sf.group_en, atlas_sf.kingdom, atlas_sf.phylum, atlas_sf.class, atlas_sf.order, atlas_sf.family, atlas_sf.genus, qc_pix.ID AS pix_id, qc_pix.geom AS geometry, to_binary(qc_pix.geom) AS geometry_text FROM atlas_sf, qc_pix WHERE atlas_sf.year_obs == " + str(year) + " AND ST_Intersects(atlas_sf.geometry, qc_pix.geom)"
    print(qu)
    d.sql(qu)
    # d.sql("SHOW atlas_pix")
    # d.sql("SELECT * FROM atlas_pix")


    # Write the new db
    # ici passe par un csv puis un parquet avec GDAL car incapable de recuperer les geometries au format blob quand on cree un parquet avec SQL
    # qu2 = "COPY (SELECT * FROM atlas_pix) TO 'atlas_pix_parquet/atlas_pix_" + str(year) + ".parquet' (FORMAT 'parquet')"
    qu2 = "COPY (SELECT * FROM atlas_pix) TO '/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_pix_parquet/atlas_pix_" + str(year) + ".csv' CSV HEADER DELIMITER ','"
    d.sql(qu2)
    d.sql("DROP TABLE atlas_pix")

    req = "ogr2ogr -f Parquet -s_srs EPSG:4326 -t_srs EPSG:4326 /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_pix_parquet/atlas_pix_" + str(year) + ".parquet /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_pix_parquet/atlas_pix_" + str(year) + ".csv"
    os.system(req)
    os.system("rm /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_pix_parquet/*.csv")


# ------------------------------------------------------------- #
# Zone de verification
# ------------------------------------------------------------- #

# Zone de test pour writer un file en geoparquet
year=2015
qu = "CREATE TABLE atlas_pix AS SELECT * FROM atlas_sf, qc_pix WHERE atlas_sf.year_obs == " + str(year) + " AND ST_Intersects(atlas_sf.geometry, qc_pix.geom)"
print(qu)
d.sql(qu)

d.sql("SELECT * FROM atlas_pix")

# Etapes de Guillaume: requete sql sur db, write in csv, converti le csv en parquet avec GDAL
# cf script geospatial-docker/atlas2pmtiles.sh 
d.sql("COPY (SELECT * FROM atlas_pix) TO 'atlas2015_test.csv' CSV HEADER DELIMITER ','")

#GDAL
ogr2ogr -f Parquet -s_srs EPSG:4326 -t_srs EPSG:4326 atlas2015_test.parquet atlas2015_test.csv

#--> voir avec Guillaume
# pb est qu'il manque REQUIRED METADATA pour valider le parquet en geoparquet : A GeoParquet file MUST include a geo key in the Parquet metadata (see FileMetaData::key_value_metadata). The value of this key MUST be a JSON-encoded UTF-8 string representing the file and column metadata that validates against the GeoParquet metadata schema. The file and column metadata fields are described below (cf https://geoparquet.org/releases/v1.0.0/ & https://github.com/geopandas/geopandas/discussions/3158) 




#### Extraction of the spe richness and spe list for each pix for each year for each taxa group ####
import geopandas as geopd
import pandas as pd
import matplotlib.pyplot as plt
from rasterio.plot import show
import os
import rasterio as rio
from rasterio import features
from rasterio.enums import MergeAlg
from rasterio.mask import mask
from numpy import int16  

init_path="/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/"

# qc = geopd.read_file(init_path + "qc_polygons/QUEBEC_unique_poly.gpkg").to_crs(4326)
# qc_json = qc.to_json()
# qc.to_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/qc_poly.json", driver="GeoJSON")

qc_pix=geopd.read_file(init_path + "qc_polygons/qc_grid_1x1km_finale_latlon.gpkg")
qc_pix=qc_pix.rename(columns={'ID':'pix_id'})

qc.plot()
plt.show()

global_spe_rich=geopd.GeoDataFrame(columns=["pix_id", "spe_rich", "spe_list", "group_en", "year"])

taxa_group = ["Unknown", "Amphibians", "Birds", "Other taxons", "Reptiles", "Fungi", "Angiosperms", "Vascular cryptogam", "Mammals", "Algae", "Other invertebrates", "Bryophytes", "Other plants", "Fish", "Arthropods", "Conifers", "Tunicates"]

for year in range(1900, 2024):
# for year in range(2020, 2024):
    print(year)
    # for taxa in taxa_group:
    #     print(year, taxa)
    data_path=init_path + "atlas_pix_parquet/atlas_pix_" + str(year) + ".parquet"
    data=pd.read_parquet(data_path)
    # uniq_pix=data.drop_duplicates("geometry")[["pix_id", "geometry"]]
    

    spe_list=data.groupby(["pix_id", "group_en"])["valid_scientific_name"].unique()
    spe_list=spe_list.reset_index() # split le multi-index in two columns
    spe_rich=data.groupby(["pix_id", "group_en"])["valid_scientific_name"].nunique()
    info_pix=pd.DataFrame({'pix_id':spe_list.pix_id, 'spe_rich':spe_rich.values, 'spe_list':spe_list.valid_scientific_name, 'group_en':spe_list.group_en, 'year':year})
    info_pix['pix_id']=pd.to_numeric(info_pix['pix_id'])

    # update GEOdataframe with species richness and list of species
    # -------------------------------------------------------------
    global_spe_rich=pd.concat([global_spe_rich, info_pix])

# create a geopackage per year per taxo
for year in range(1900, 1901):
    taxa_pres = global_spe_rich.group_en[global_spe_rich['year'] == year].unique()
    print(taxa_pres)
    for group in taxa_pres:
        # df = global_spe_rich[(global_spe_rich['year'] == year) & (global_spe_rich['group_en'] == group)]
        print(group)
    

    # rasterization for each year for each taxonomic group
    # ---------------------------------------------------

    # 1 - left join to get information for all pixels per taxo group
    # --------------------------------------------------------------
    # data2 = info_pix.drop(columns=["spe_list"])

    # for group in taxa_group:
            
    #     data_group = data2[data2["group_en"] == group]
    #     data_group_join = qc_pix.merge(data_group, how='left', on='pix_id')
    #     data_group_join['spe_rich'] = data_group_join['spe_rich'].fillna(0)
    #     data_group_join['year'] = data_group_join['year'].fillna(year)
        
    #     data_group_join.plot(column = "spe_rich")
    #     plt.show()
# test avec utilisation geopackage
        # data_group_join.to_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/qc_pix_rs_all_taxa.gpkg", driver = "GPKG")







#### ===> reprendre ici pour la rasterisation
        # shapes = ((geom,value) for geom, value in zip(data_group_join.geometry, data_group_join["spe_rich"]))

        # raster = rasterio.features.rasterize(shapes = shapes, out_shape = )
        # data_group_join.spe_rich.min()
        # data_group_join.spe_rich.max()
        # data_group_join.spe_rich.describe()

        # 2 - conversion to geodataframe
        # ------------------------------
        group_sf=geopd.GeoDataFrame(data_group_join, geometry = "geometry", crs = "EPSG:4326")

        # 3 - from geodf to geojson
        # -------------------------
        # geoj=data2.to_json()
        # group_sf.to_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/geojson_test.json", driver="GeoJSON")

        # 4 - from geojson to raster with rasterio
        # ----------------------------------------
        # os.system("rio rasterize /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/test.tif --res 0.0167 < /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/geojson_test.json")

        # 3 - visualisation
        # raster1=rio.open("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/test.tif")
        # raster1.meta
        # plt.imshow(raster1.read(1))
        # plt.show()

    # final_sf.plot(column='spe_rich')
    # plt.show()

    print("----------> Year " + str(year) + " DONE !")





# Visual explo
fig, (ax1, ax2) = plt.subplots(1,2)
show((raster1, 1), ax=ax1)
# qc.plot(ax=ax1, facecolor="none", edgecolor="grey", linewidth = 0.25)
show((raster1_crop, 1), ax=ax2)
# qc.plot(ax=ax2, facecolor="none", edgecolor="grey", linewidth = 0.25)
plt.show()

# Pour calcul sur raster, necessaire d'enlever les nan
import numpy as np
rast2=raster1[~np.isnan(raster1)]
rast2.sum()
rast2.max()


import matplotlib.pyplot as plt

plt.imshow(raster1_crop, cmap='viridis')
plt.colorbar(label='Data Values')
plt.title('Raster with Nodata Values')
plt.show() 

# Exploration
global_spe_rich.dtypes
global_spe_rich.spe_rich.describe()
max(global_spe_rich.spe_rich)

# In order to create a geoparquet
# see https://geoparquet.org/releases/v1.1.0/schema.json
# and https://geoparquet.org/releases/v1.1.0/
global_spe_rich.index
global_spe_rich.to_parquet(path="/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_pix_parquet/EBV_rich_spe_raw_data.parquet")
                        #    , index=False, compression="None", geometry_encoding="WKB", write_covering_bbox=False, schema_version=None)
# *** write_covering_bbox: Writes the bounding box column for each row entry with column name ‘bbox’. Writing a bbox column can be computationally expensive, but allows you to specify a bbox in : func:read_parquet for filtered reading. Note: this bbox column is part of the newer GeoParquet 1.1 specification and should be considered as experimental. While writing the column is backwards compatible, using it for filtering may not be supported by all readers.
# *** aactual version of geoparquet on my laptop is 0.0.3, not accepted for the argument schema_version ({‘0.1.0’, ‘0.4.0’, ‘1.0.0’, ‘1.1.0’, None})
# pip freeze pour voir les versions des paquets installés

# test parquet --> GOOD !
data_path="/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/atlas_pix_parquet/EBV_rich_spe_raw_data.parquet"
rich=pd.read_parquet(data_path)
rich["geometry"] = geopd.GeoSeries.from_wkb(rich["geometry"])
rich_sf=geopd.GeoDataFrame(rich, geometry = "geometry", crs = "EPSG:4326")
rich_sf





rs = geopd.read_file("/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test/geojson_test.json")