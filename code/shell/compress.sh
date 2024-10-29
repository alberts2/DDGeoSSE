#!/usr/bin/env sh

# compress the output files in the remote server

zip param_label.zip sim.*.labels.csv #labeled parameters 
zip anagen_changes.zip anagenetic_changes_tree_*.csv #anagenesis history over time
find . -name "species_count_region_tree_*.csv" -print0 | xargs -0 zip species_count_reg.zip #species count in regions over time
find . -name "species_count_range_tree_*.csv" -print0 | xargs -0 zip species_count_range.zip #species count in ranges over time
find . -name "cts_species_count_region_tree_*.csv" -print0 | xargs -0 zip cts_species_count_reg.zip #species count in regions over cts time
find . -name "cts_species_count_range_tree_*.csv" -print0 | xargs -0 zip cts_species_count_range.zip #species count in ranges over cts time
find . -name "sim.*.dat.csv" -print0 | xargs -0 zip tip_data.zip #character data at tips
find . -name "sim.*.tre" -print0 | xargs -0 zip complete_tree.zip #complete trees 
find . -name "cladogenetic event history_*.csv" -print0 | xargs -0 zip cladogen_history.zip #cladogenesis history over time
