# This script will simulate phylogenetic datasets under a 
# density-dependendent GeoSSE model.

#####################
# Load julia packages
#####################
using Combinatorics # Package for combinatorics/permutations
using OrderedCollections # Package for ordering dictionary 
using DataStructures # Package for preserving dictionary order by its insertion order 
using Distributions # Package for using probability distributions 


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
num_ranges  = 2^num_regions - 1; # Number of possible ranges given num_regions
num_species_per_range  = fill(0, num_ranges); # Empty vector to store number of species in each range
# num_species_per_region = fill(0, num_regions); # Empty vector to store number of species in each region
num_species_per_region = [1,0,0] # just for testing (remove later)
# P.S. the entries for num_species_per_region and num_species_per_range will be updated at each time step

#################################################
# Create dictionary for region occupancy status (present/absent) and range labels
#################################################
# e.g. (1,0,0) means a species present only in region A, so it has range {A}.
# ranges_inv = {
#     (1,0,0):1,
#     (0,1,0):2,
#     (0,0,2):3,
#     (1,1,0):4,
#     (1,0,1):5,
#     (0,1,1):6,
#     (1,1,1):7
# }
dict_values = zeros(Int64,0); # Intialize empty vector for the values field of the dictionary
for i in 1:num_ranges
    push!(dict_values,i); # Populate the values field 
end

dict_keys = Vector{Int}[]; # Create an empty vector of vectors for the keys field of the dictionary
for j in 1:num_regions
  zero_vec = zeros(Int64,num_regions); # Create a zero vector of length num_regions
  zero_vec[1:j] = fill(1,j); # Populating the vector with 1's for indicating presence/absence of species
  range_comb = unique(permutations(zero_vec)); # different permutations of (1,...) for range size <= num_regions 
  for k in 1:length(range_comb)
    push!(dict_keys,range_comb[k]); # populating keys field of the dictionary
  end
end

ranges_inv = Dict(dict_keys .=> dict_values) # Append the keys and values to the dictionary 
ranges_inv = sort(ranges_inv; byvalue=true) # sort the dictionary by its values field
# ranges = [ [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [1,1,0], [1,1,1] ]


# Define event space as dictionaries (for checking if we have a valid cladogenesis/anagenesis during simulation)
################################################
# Dictionary for within-region speciation (evt_w)
################################################
# key is ancestral range, value is list of possible child ranges that differ from ancestral range
# for example:
#     (1,1,0):[(1,0,0),(0,1,0)],
# valid events:
#     AB->AB,A (budding cladogenesis)
#     AB->A,AB (budding cladogenesis)
#     AB->AB,B (budding cladogenesis)
#     AB->B,AB (budding cladogenesis)
# evt_w = {
#     (1,0,0):[(1,0,0)],
#     (0,1,0):[(0,1,0)],
#     (0,0,1):[(0,0,1)],
#     (1,1,0):[(1,0,0),(0,1,0)],
#     (1,0,1):[(1,0,0),(0,0,1)],
#     (0,1,1):[(0,1,0),(0,0,1)],
#     (1,1,1):[(1,0,0),(0,1,0),(0,0,1)]
# }

evt_w       = OrderedDict(); # Create an empty dictionary for within-region speciation 
dict_keys_w = copy(dict_keys); # keys field of the dictionary evt_w

for k in 1:length(dict_keys_w)
    if count(==(1),dict_keys_w[k])==1 # if the ancestral range is of size 1 (endemic)
        evt_w[dict_keys_w[k]] = dict_keys_w[k]; # assign possible range for children = ancestral range 
    elseif count(==(1),dict_keys_w[k])>1 # if the ancestral range is of size > 1 (widespread)
        idx = findall(dict_keys_w[k].== 1 ); # find locations of 1's in ancestral species range k
        child_range_w = Vector{Int}[]; # initialize the vector of vectors to store all possible child ranges from ancestral of range k
        for l in 1:length(idx) # loop over each 1 in ancestral species range k for within-region speciation in that location
            child_range_idx = zeros(Int64,num_regions)
            child_range_idx[idx[l]] = 1 # indicates that speciation happens in location idx[l] from the parent species with range k
            push!(child_range_w,child_range_idx) # append to possible range for child species from ancenstral species with range k
        end
        evt_w[dict_keys_w[k]] = child_range_w; # create dictionary for evt_w
    end
end

################################################
# Dictionary for extinction (evt_e)
################################################
# key is range before extinction event, value is all possible ranges following an extinction. 
# [(0,0,0)] indicates species extinction
#
# evt_e = {
#     (1,0,0):[(0,0,0)],
#     (0,1,0):[(0,0,0)],
#     (0,0,1):[(0,0,0)],
#     (1,1,0):[(1,0,0),(0,1,0)],
#     (1,0,1):[(1,0,0),(0,0,1)],
#     (0,1,1):[(0,1,0),(0,0,1)],
#     (1,1,1):[(1,1,0),(1,0,1),(0,1,1)]
# }
evt_e       = OrderedDict(); # Create an empty dictionary for extinction
dict_keys_e = copy(dict_keys); # keys field of the dictionary evt_e

for k in 1:length(dict_keys_e)
    if count(==(1),dict_keys_e[k])==1 # if the range is of size 1 (endemic)
        evt_e[dict_keys_e[k]] = zeros(Int64,num_regions); # species extinction 
    elseif count(==(1),dict_keys_e[k])>1 # if the ancestral range is of size > 1 (widespread)
        idx = findall(dict_keys_e[k].== 1 ); # find locations of 1's in species range k before extinction
        new_range_e = Vector{Int}[]; # initialize the vector of vectors to store all possible new ranges after extinction
        for l in 1:length(idx) # loop each location index with a species in range k
            new_range_idx = copy(dict_keys_e[k])
            new_range_idx[idx[l]] = 0 # indicates that extinction occurs in location idx[l] in species with range k
            push!(new_range_e,new_range_idx) # append to possible range new range after extinction for species with range k
        end
        evt_e[dict_keys_e[k]] = new_range_e; # create dictionary for evt_e
    end
end

################################################
# Dictionary for dispersal (evt_d)
################################################
# # key is range before dispersal, value indicates all possible ranges following a dispersal event 
# evt_d = {
#     (1,0,0):[(1,1,0),(1,0,1)],
#     (0,1,0):[(1,1,0),(0,1,1)],
#     (0,0,1):[(0,1,1),(0,1,1)],
#     (1,1,0):[(1,1,1)],
#     (1,0,1):[(1,1,1)],
#     (0,1,1):[(1,1,1)],
#     (1,1,1):[(,)]
# }

evt_d       = OrderedDict(); # Create an empty dictionary for dispersal
dict_keys_d = copy(dict_keys); # keys field of the dictionary evt_d

for k in 1:length(dict_keys_d)
    if count(==(1),dict_keys_d[k])==num_regions # if the species is already widespread in all regions 
        evt_d[dict_keys_d[k]] = [] # can't disperse
    elseif count(==(1),dict_keys_d[k])!=num_regions # if the range is still < num_regions, can disperse
        idx = findall(dict_keys_d[k].== 0 ); # find locations of all 0's in species with range k before dispersal
        new_range_d = Vector{Int}[]; # intialize vector of vectors to store all possible new ranges following a dispersal of species with range k
        for l in 1:length(idx) # loop each location index with a species in range k
            new_range_idx = copy(dict_keys_d[k])
            new_range_idx[idx[l]] = 1 # indicates that species dispersal to that location idx[l]
            push!(new_range_d,new_range_idx) # append to possible range new range after extinction for species with range k
        end
        evt_d[dict_keys_d[k]] = new_range_d; # create dictionary for evt_d
    end
end


################################################
# Dictionary for between-region speciation (evt_b)
################################################
# key indicate ancestral species range, value indicates all possible child ranges for one of the two children 
# following the speciation event. List of all possible ranges of the other child can be computed using:
# ancestral range - child range. E.g. ancestral = (1,1,0), child = (0,1,0). Then the other child is: 
# (1,1,0) - (0,1,0) = (1,0,0)

# evt_b = {
#     (1,0,0):[(,)],
#     (0,1,0):[(,)],
#     (0,0,1):[(,)],
#     (1,1,0):[(1,0,0)],
#     (1,0,1):[(1,0,0)],
#     (0,1,1):[(0,1,0)],
#     (1,1,1):[(1,0,0),(0,1,0),(0,0,1)]
# }

evt_b       = OrderedDict(); # Create an empty dictionary for between-region speciation
dict_keys_b = copy(dict_keys); # keys field of the dictionary evt_b

for k in 1:length(dict_keys_b) # Looping over each ancestral species range in keys field 
    if count(==(1),dict_keys_b[k])==1 # if the ancestral species is endemic, can't speciate 
        evt_b[dict_keys_b[k]] = [] # then assign an empty vector for child range
    # the special case when ancestral range of size 2, since we can only have 1|1 split. So, no need to iterate all permutations of (1,0,0,...) 
    elseif count(==(1),dict_keys_b[k])==2 
        idx_anc = findall(dict_keys_b[k].== 1 ); # find locations of 1's in ancestral range k with size 2
        child_range_b = zeros(Int64,num_regions); # Initialize vector of zeros 
        child_range_b[idx_anc[1]] = 1; # Assign 1
        evt_b[dict_keys_b[k]] = child_range_b;
    else # if the ancestral species is widespread, can speciate
        anc_range_size = count(==(1),dict_keys_b[k]) # ancestral range size for index k in key field
        # looping for each different size of child range after speciation.
        # Divide by 2 because the range size of the other child can be derived from the range of the other one. 
        idx_anc = findall(dict_keys_b[k].== 1 ); # find locations of 1's in ancestral range k
        #
        child_range_b = Vector{Int}[]; # initialize vector of vectors to store all child ranges from ancestral range k 
        for l in 1:floor.(Int64,anc_range_size/2) # loop over each range split size for the child range e.g. for k = 4, we can have split of size 1&2 (not including its compelement of size 3)
            child_init = zeros(Int64,num_regions); # Initialize vector of zeros 
            child_init[1:l] .= 1; # Assign 1's equal to the range split size 
            child_comb = unique(permutations(child_init)); # permute all possibilities of combinations of 1's for that split size 
            for m in 1:length(child_comb) # loop over this permutation 
                idx_child = findall(child_comb[m].== 1 );
                test_idx = [x in idx_anc for x = idx_child] # Check if that split of size l is valid (i.e. if its contained in ancestral range)
                if all(test_idx) == true
                    push!(child_range_b,child_comb[m]) # If yes, then we accept that split of size l as a valid split
                end
            end
        end
        evt_b[dict_keys_b[k]] = child_range_b; # append keys and value fields to dictionary 
    end
end 
# println(evt_b)
# println(length(evt_b))
# xxx

###########################################
# # Define rates for event space in each model variant 
###########################################

#########################
# Rate function for within-region speciation 
#########################

# Rate function for density-independent within-region speciation

function rate_w(base,region_ind,rate_scalar)
    ## Error checking
    @assert base > 0 "base rate must be positive" # base rate must be positive 
    @assert region_ind in dict_values "region index must be <= num_ranges" # region ind must be element of dict_values (region labels)
    n_i   = num_species_per_region[region_ind];  @assert n_i > 0 "number of species in region_ind must be > 0"
    ## Define rate parameters 
    rho_w = base # base rate 
    w_d   = rate_scalar # rate scalar for within-region speciation
    m_i_w = rand(Uniform(0.01,1)) # draw density-independent rate randomly [0.01,1]
    ## Compute the overall rates
    m_d_w = (log(n_i)+1)^(w_d)
    ans   = base*m_i_w*m_d_w
    return(ans)
end
# test = rate_w(0.1,1,0.1)
# println(test)
# xxx


# # simulate
# do_simulation = true
# start_time = 0.0
# start_state = runif(1, 1, num_ranges)


# while (do_simulation) {
    
#     # draw time to next event
    
#     # draw next event type
    
#     # update current simulation state
    
#     # check if we stop simulation
#     if stop_condition == "time" and curr_time > max_time {
#         do_simulation = false
#     } else if stop_condition == "taxa" and curr_taxa > max_taxa {
#         do_simulation = false
#     }
    
#     # stop simulation if condition met
#     if (!do_simulation) {
#         break
#     }
    
# }

# # output

# # write tree to file

# # write ranges to file

# # write parameters to file

# # write summary statistics to file