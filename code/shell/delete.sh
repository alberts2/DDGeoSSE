#!/bin/bash

# remove output files from the remote server 

for i in $(seq 1 50000); do
    rm sim.$i.labels.csv   
    rm anagenetic_changes_tree_$i.csv 
    rm species_count_region_tree_$i.csv
    rm species_count_range_tree_$i.csv 
    rm cts_species_count_region_tree_$i.csv
    rm cts_species_count_range_tree_$i.csv 
    rm sim.$i.dat.csv
    rm sim.$i.tre
    rm "cladogenetic event history_$i.csv"
done

rm anagen_changes.zip
rm complete_tree.zip 
rm cts_species_count_range.zip
rm cts_species_count_reg.zip
rm param_label.zip
rm species_count_range.zip
rm species_count_reg.zip
rm tip_data.zip
rm cladogen_history.zip