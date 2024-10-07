#!/bin/bash

for i in $(seq 1 50000); do
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/param_label_tree_$i.labels.csv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/anagenetic_changes_tree_$i.csv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/species_count_region_tree_$i.csv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/species_count_range_tree_$i.csv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/cts_species_count_region_tree_$i.csv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/cts_species_count_range_tree_$i.csv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/tip_data_tree_$i.dat.csv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
    mv /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/code/jl/complete_tree_$i.tre /Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/
done