# This shell script simulates batches of trees 
# sh sim_batch_phyddle.sh 2 2 simulate sim true 1000

#!/usr/bin/env sh

# Number of batches
#batches=500
batches=$1
NUMSIM=$2

out_path=$3
out_prefix=$4
resimulate=$5
timebin=$6


# Number of parallel processes (cores) per batch
cores=2

# seq 1 $batches | xargs -n 1 -P $cores -I {} sh -c './sim_one.sh {}' 
seq 1 $batches | xargs -P $cores -I {} sh -c "./sim_one_phyddle.sh {} ${NUMSIM} ${out_path} ${out_prefix} ${resimulate} ${timebin}"


# Wait for all background processes to finish
wait
echo "All batches are complete."


