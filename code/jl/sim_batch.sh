# This shell script script simulate batches of batch_size simulation replicates in parallel under sim_one_func.jl

#!/usr/bin/env sh

# Number of batches
batches=10

# Number of simulation replicates per batch
batch_size=1000

# Number of parallel processes (cores) per batch
cores=10

# Julia arguments
out_path=simulate
out_prefix=sim

seq 1 $batches | xargs -P $cores -I {} sh -c 'batch_num={}
    start_idx=$(( (batch_num - 1) * '"${batch_size}"' + 1 ))
    julia --quiet sim_one_func.jl '"${out_path}"' '"${out_prefix}"' "${start_idx}" '"${batch_size}"' > /dev/null 2>&1
'
# Wait for all background processes to finish
wait
echo "All batches are complete."

