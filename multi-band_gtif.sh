#!/bin/bash
for year in {1992..2017..1}
do
    export YEARBLOC=${year} #export the year to be used in Rscript

    s5cmd cp 's3://bq-io/acer/oiseaux-nicheurs-qc/*range_'${year}'*' /home/claire/desktop/data/
    gdal_merge.py -separate -o /home/claire/desktop/data/multiband/${year}_range_multiband.tif /home/claire/desktop/data/*.tif

    Rscript  /home/claire/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/multi-band_gtif.r

    rm /home/claire/desktop/data/*.tif


done