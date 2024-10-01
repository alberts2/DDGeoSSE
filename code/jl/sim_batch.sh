# This shell script simulates batches of trees 

#!/usr/bin/env sh

# Number of batches
#batches=500
batches=2

# Number of parallel processes (cores) per batch
cores=2

# seq 1 $batches | xargs -n 1 -P $cores -I {} sh -c './sim_one.sh {}' 
seq 1 $batches | xargs -P $cores -I {} sh -c './sim_one.sh {}' 

# Wait for all background processes to finish
wait
echo "All batches are complete."


