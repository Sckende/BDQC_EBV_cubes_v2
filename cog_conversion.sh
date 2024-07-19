#!/bin/bash

# cd /home/local/USHERBROOKE/juhc3201/Desktop/rs_inla

# for i in {1992..2017..1}
# do
#     gdaladdo -r average path 2 4 8 16 
#     gdal_translate rs_inla_${i}.tif rs_inla_${i}_cog.tif -co COMPRESS=LZW -co TILED=YES -co COPY_SRC_OVERVIEWS=YES 
# done


INFILE=/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/INLA_file_name.txt

# cat ${INFILE}
file_path=/home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/data_test

# IFS=$'\n'

for LINE in $(cat ${INFILE})
do
    # echo ${LINE}
    # echo ${file_path}/cog/${LINE/.tif/_cog.tif}

    gdaladdo -r average ${file_path}/${LINE} 2 4 8 16 
    gdal_translate ${file_path}/${LINE} ${file_path}/cog/${LINE/.tif/_cog.tif} -co COMPRESS=LZW -co TILED=YES -co COPY_SRC_OVERVIEWS=YES 
done