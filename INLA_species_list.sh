#!/bin/bash

# Pour obtenir la liste des especes d'oiseaux etudie par Vincent Bellavance avec INLA

# Step 1 
# creation d'un ficher pour récupérer la liste des modeles à partir du clooud

s5cmd ls 's3://bq-io/acer/TdeB_benchmark_SDM/oiseaux-nicheurs-qc/*range_2017*' > /home/claire/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/INLA_species_list_init.txt

# step 2 - traitement pour obtenir uniquement la liste des especes
awk '{print $4}' INLA_species_list.txt | awk -F'_' '{print $1"_"$2}' > /home/claire/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/INLA_species_list_final.txt
