# This script will simulate phylogenetic datasets under a 
# density-dependendent GeoSSE model.


# input
base_dir = "~/projects/density_geosse"

# settings
num_regions = 3
min_taxa = 5
max_taxa = 100
max_time = 10.0
start_state = 0    # or draw randomly
stop_condition = "time"

# set up model
# define state space
num_regions = 3
num_ranges = 2^3 - 1
num_species_per_range = repeat(0, num_ranges)
num_species_per_region = repeat(0, num_states)
# ranges for A, B, C, AB, AC, BC, ABC
ranges = [ [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [1,1,0], [1,1,1] ]
ranges_inv = {
    (1,0,0):1,
    (0,1,0):2,
    (0,0,2):3,
    (1,1,0):4,
    (1,0,1):5,
    (0,1,1):6,
    (1,1,1):7
}

# define event space

# Within-region speciation
# key is ancestral range, value is list of possible child ranges that differ from ancestral range
# for example:
#     (1,1,0):[(1,0,0),(0,1,0)],
# valid events:
#     AB->AB,A
#     AB->A,AB
#     AB->AB,B
#     AB->B,AB
evt_w = {
    (1,0,0):[(1,0,0)],
    (0,1,0):[(0,1,0)],
    (0,0,1):[(0,0,1)],
    (1,1,0):[(1,0,0),(0,1,0)],
    (1,0,1):[(1,0,0),(0,0,1)],
    (0,1,1):[(0,1,0),(0,0,1)],
    (1,1,1):[(1,0,0),(0,1,0),(0,0,1)]
}

# Extinction
# ...
evt_e = {
    (1,0,0):[(0,0,0)],
    (0,1,0):[(0,0,0)],
    (0,0,1):[(0,0,0)],
    (1,1,0):[(1,0,0),(0,1,0)],
    (1,0,1):[(1,0,0),(0,0,1)],
    (0,1,1):[(0,1,0),(0,0,1)],
    (1,1,1):[(1,1,0),(1,0,1),(0,1,1)]
}

# Dispersal
# ...
evt_d = {
    (1,0,0):[(1,1,0),(1,0,1)],
    (0,1,0):[(1,1,0),(0,1,1)],
    (0,0,1):[(0,1,1),(0,1,1)],
    (1,1,0):[(1,1,1)],
    (1,0,1):[(1,1,1)],
    (0,1,1):[(1,1,1)],
    (1,1,1):[(,)]
}

# Between-region speciation
evt_b = {
    (1,0,0):[(,)],
    (0,1,0):[(,)],
    (0,0,1):[(,)],
    (1,1,0):[(1,0,0)],
    (1,0,1):[(1,0,0)],
    (0,1,1):[(0,1,0)],
    (1,1,1):[(1,1,0),(1,0,1),(0,1,1)]
}

# define rates for event space

# simulate
do_simulation = true
start_time = 0.0
start_state = runif(1, 1, num_ranges)


while (do_simulation) {
    
    # draw time to next event
    
    # draw next event type
    
    # keep track of phylogeny and parent-child relationships as it evolves
    
    # update current simulation state
    
    # check if we stop simulation
    if stop_condition == "time" and curr_time > max_time {
        do_simulation = false
    } else if stop_condition == "taxa" and curr_taxa > max_taxa {
        do_simulation = false
    }
    
    # stop simulation if condition met
    if (!do_simulation) {
        break
    }
    
}

# output

# write tree to file

# write ranges to file

# write parameters to file

# write summary statistics to file