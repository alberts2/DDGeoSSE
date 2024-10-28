#!/usr/bin/env sh

# compress the output files in the remote server

zip param_label.zip param_label_tree_*.labels.csv
zip anagen_changes.zip anagenetic_changes_tree_*.csv
find . -name "species_count_region_tree_*.csv" -print0 | xargs -0 zip species_count_reg.zip
find . -name "species_count_range_tree_*.csv" -print0 | xargs -0 zip species_count_range.zip
find . -name "cts_species_count_region_tree_*.csv" -print0 | xargs -0 zip cts_species_count_reg.zip
find . -name "cts_species_count_range_tree_*.csv" -print0 | xargs -0 zip cts_species_count_range.zip
find . -name "tip_data_tree_*.dat.csv" -print0 | xargs -0 zip tip_data.zip
find . -name "complete_tree_*.tre" -print0 | xargs -0 zip complete_tree.zip
find . -name "cladogenetic event history_*.csv" -print0 | xargs -0 zip cladogen_history.zip
