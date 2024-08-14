#!/bin/bash
mkdir /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/dl_acer/
for year in {1992..2017..1}
do
    export YEARBLOC=${year} #export the year to be used in Rscript

    s5cmd cp 's3://bq-io/acer/oiseaux-nicheurs-qc/*range_'${year}'*' /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/dl_acer/

    gdal_merge.py -separate -o /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/inla_multiband/${year}_range_multiband.tif /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/dl_acer/*.tif

    Rscript  /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/multi-band_gtif.r

    rm /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/dl_acer/*.tif


done

rm -rf /home/local/USHERBROOKE/juhc3201/BDQC-GEOBON/data/QUEBEC_in_a_cube/Richesse_spe_version_2/dl_acer/