# This shell script script simulate NUMSIM trees under batch BATCH_NO

#!/usr/bin/env sh
BATCH_NO=$1

#NUMSIM=100
NUMSIM=2

timebin=1000

resimulate="true"


# Run 1 batch of simulations
for IDX in $(seq $(((BATCH_NO-1)*NUMSIM+1)) $((BATCH_NO*NUMSIM)));do
    # julia --quiet sim_test_v10.jl ${IDX} ${resimulate} ${timebin} > /dev/null 2>&1
    julia --quiet sim_test_v11.jl ${IDX} ${resimulate} ${timebin} > /dev/null 2>&1 
    # echo "batch ${BATCH_NO} tree ${IDX} complete"
done

