#!/bin/bash

# remove output files from the remote server 

for i in $(seq 1 50000); do
    rm param_label_tree_$i.labels.csv   
    rm anagenetic_changes_tree_$i.csv 
    rm species_count_region_tree_$i.csv
    rm species_count_range_tree_$i.csv 
    rm cts_species_count_region_tree_$i.csv
    rm cts_species_count_range_tree_$i.csv 
    rm tip_data_tree_$i.dat.csv
    rm complete_tree_$i.tre
done

rm anagen_changes.zip
rm complete_tree.zip 
rm cts_species_count_range.zip
rm cts_species_count_reg.zip
rm param_label.zip
rm species_count_range.zip
rm species_count_reg.zip
rm tip_data.zip