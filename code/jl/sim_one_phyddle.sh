# This shell script script simulate NUMSIM trees under batch BATCH_NO

#!/usr/bin/env sh
BATCH_NO=$1

NUMSIM=$2

out_path=$3
out_prefix=$4
resimulate=$5
timebin=$6

echo "BATCH_NO: $BATCH_NO"
echo "NUMSIM: $NUMSIM"
echo "out_path: $out_path"
echo "out_prefix: $out_prefix"
echo "resimulate: $resimulate"
echo "timebin: $timebin"

# Run 1 batch of simulations
for IDX in $(seq $(((BATCH_NO-1)*NUMSIM+1)) $((BATCH_NO*NUMSIM)));do
    # julia --quiet sim_test_v11.jl ${IDX} ${resimulate} ${timebin} > /dev/null 2>&1 
    # julia --quiet sim_test_v12.jl ${IDX} ${resimulate} ${timebin} > /dev/null 2>&1
    julia --quiet sim_one.jl ${out_path} ${out_prefix} ${IDX} ${resimulate} ${timebin} > /dev/null 2>&1
    # echo "batch ${BATCH_NO} tree ${IDX} complete"
done

