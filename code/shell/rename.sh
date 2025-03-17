#!/bin/bash

# Initialize counter for new file indices
new_idx=350001
new_idx_2=350001
new_idx_3=350001
new_idx_4=350001

# Loop over files sim.1.tre to sim.50000.tre
find . -maxdepth 1 -name "sim.*.tre" ! -name "*.downsampled.tre" ! -name "*.extant.tre" | sort -V | while read -r old_file; do
    cp "$old_file" "../../combined_dat/simulate/sim.${new_idx}.tre"
    ((new_idx++))
done

# Loop over files sim.1.labels.csv to sim.50000.labels.csv
find . -maxdepth 1 -name "sim.*.labels.csv" | sort -V | while read -r old_file; do
    cp "$old_file" "../../combined_dat/simulate/sim.${new_idx_2}.labels.csv"
    ((new_idx_2++))
done

# Loop over files sim.*.dat.csv (excluding .extant.dat.csv)
find . -maxdepth 1 -name "sim.*.dat.csv" ! -name "*.extant.dat.csv" | sort -V | while read -r old_file; do
    cp "$old_file" "../../combined_dat/simulate/sim.${new_idx_3}.dat.csv"
    ((new_idx_3++))
done

#Loop over files sim.*.extant.dat.csv (only those ending with .extant.dat.csv)
find . -maxdepth 1 -name "sim.*.extant.dat.csv" | sort -V | while read -r old_file; do
    cp "$old_file" "../../combined_dat/simulate/sim.${new_idx_4}.extant.dat.csv"
    ((new_idx_4++))
done