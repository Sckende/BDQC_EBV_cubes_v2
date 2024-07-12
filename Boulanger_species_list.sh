#!/bin/bash

# Pour obtenir la liste des especes d'arbres etudiees par Boulanger

# Step 1 
# creation d'un ficher pour récupérer la liste des modeles à partir du clooud

s5cmd ls 's3://bq-io/acer/TdeB_benchmark_SDM/Boulanger_map_converted/aggregate_*' > /home/claire/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/Boulanger_species_list_init.txt

# step 2 - traitement pour obtenir uniquement la liste des especes
awk '{print $4}' /home/claire/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/Boulanger_species_list_init.txt | awk -F'_' '{print $2}' > /home/claire/BDQC-GEOBON/GITHUB/BDQC_EBV_cubes_v2/Boulanger_species_list_final.txt
