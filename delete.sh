#!/bin/bash

for i in $(seq 1 50000); do
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/param_label_tree_$i.labels.csv   
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/anagenetic_changes_tree_$i.csv 
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/species_count_region_tree_$i.csv
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/species_count_range_tree_$i.csv 
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/cts_species_count_region_tree_$i.csv
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/cts_species_count_range_tree_$i.csv 
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/tip_data_tree_$i.dat.csv
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/complete_tree_$i.tre
    rm /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/extant_tree_$i.tre
done