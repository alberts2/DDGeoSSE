# This script will simulate phylogenetic datasets under a 
# density-dependendent GeoSSE model.

#####################
# Load julia packages
#####################
using Combinatorics # Package for combinatorics/permutations
using OrderedCollections # Package for ordering dictionary 

#############
# File system
#############
base_dir = "~/projects/density_geosse";
out_dir  = base_dir*"/output";

################
# Intial setting
################
min_taxa       = 5; # Minimum number of extant taxa to simulate 
max_taxa       = 100; # Maximum number of extant taxa to simulate
max_time       = 10.0; # Maximum simulation time 
start_state    = 0;        # or draw randomly
stop_condition = "time" # If time is the stopping condition
#stop_condition = "numtaxa" # If the final number of extant taxa is the stopping condition

#############
# Model setup
#############
# Define state space 
num_regions = 3; # Number of regions 
num_ranges  = 2^3 - 1; # Number of possible ranges given num_regions
num_species_per_range  = fill(0, num_ranges); # Empty vector to store number of species in each range
num_species_per_region = fill(0, num_regions); # Empty vector to store number of species in each region
# Create dictionary for region occupancy status (present/absent) and range labels
# e.g. (1,0,0) means a species present only in region A, so it has range {A}.
dict_values = zeros(Int64,0); # Intialize empty vector for the values vector of the dictionary
for i in 1:num_ranges
    push!(dict_values,i); # Populate the values field 
end

dict_keys = Vector{Int}[] # Create an empty vector of vectors
for j in 1:num_regions
  zero_vec = zeros(Int64,num_regions); # Create a zero vector of length num_regions
  zero_vec[1:j] = fill(1,j); # Populating the vector with 1's for indicating presence/absence of species
  range_comb = unique(permutations(zero_vec)); # different permutations of (1,...) for range size <= num_regions 
  for k in 1:length(range_comb)
    push!(dict_keys,range_comb[k]); # populating keys field of the dictionary
  end
end

ranges_inv = Dict(dict_keys .=> dict_values)
ranges_inv = sort(ranges_inv; byvalue=true) # sort the dictionary by its values field
# ranges = [ [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [1,1,0], [1,1,1] ]
# ranges_inv = {
#     (1,0,0):1,
#     (0,1,0):2,
#     (0,0,2):3,
#     (1,1,0):4,
#     (1,0,1):5,
#     (0,1,1):6,
#     (1,1,1):7
# }

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