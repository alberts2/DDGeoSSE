# This script will simulate phylogenetic datasets under a 
# density-dependent GeoSSE model.

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
# Initial setting
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
# num_regions = 4; # Number of regions 
num_regions = 3; # Number of regions 
num_ranges  = 2^num_regions - 1; # Number of possible ranges given num_regions
num_species_per_range  = fill(0, num_ranges); # Empty vector to store number of species in each range
num_species_per_region = fill(0, num_regions); # Empty vector to store number of species in each region
# P.S. the entries for num_species_per_region and num_species_per_range will be updated at each time step
# since rates only depend on n_i at the current time, not past times. 

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
        # evt_b[dict_keys_b[k]] = child_range_b;
        evt_b[dict_keys_b[k]] = [child_range_b];
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
# println(evt_b[[1,1,1,0]])
# println([1,1,1,0]-evt_b[[1,1,1,0]][1])
# println(length(evt_b))
# xxx

###########################################
# # Define rates for event space in each model variant 
###########################################

#########################
# Rate function for within-region speciation 
#########################

function rate_w(base,region_ind,rate_scalar,model)
    ## Error checking
    @assert base > 0 "base rate must be positive" # base rate must be positive 
    @assert region_ind <= num_regions && region_ind > 0 "region index must be in {1,2,...,num_regions}" # region ind must be <= num_regions
    n_i   = num_species_per_region[region_ind] # Number of species currently in region_ind
    if n_i == 0
        ans = 0 # rate equals 0 if no species 
    else
        ## Define rate parameters 
        rho_w = base # base rate 
        w_d   = rate_scalar # rate scalar for within-region speciation
        m_i_w = rand(Uniform(0.01,1)) # draw density-independent rate randomly [0.01,1]
        ## Compute the overall rates
        if model == "log"
            m_d_w = (log(n_i)+1)^(w_d)
        elseif model == "standard"
            m_d_w = n_i^w_d
        end
        ans   = rho_w*m_i_w*m_d_w
    end
    return(ans)
end
# num_species_per_region = [1,0,0] # just for testing (remove later)
# test = rate_w(0.1,1,0.1,"log")
# println(test)
# xxx

#########################
# Rate function for extinction
#########################

function rate_e(base,region_ind,rate_scalar,model)
    ## Error checking
    @assert base > 0 "base rate must be positive" # base rate must be positive 
    @assert region_ind <= num_regions && region_ind > 0 "region index must be in {1,2,...,num_regions}" # region ind must be <= num_regions
    n_i   = num_species_per_region[region_ind] # Number of species currently in region_ind
    if n_i == 0
        ans = 0 # rate equals 0 if no species 
    else # if there is species 
        ## Define rate parameters 
        rho_e = base # base rate 
        e_d   = rate_scalar # rate scalar for extinction
        m_i_e = rand(Uniform(0.01,1)) # draw density-independent rate randomly [0.01,1]
        ## Compute the overall rates
        if model == "log"
            m_d_e = (log(n_i)+1)^(e_d)
        elseif model == "standard"
            m_d_e = n_i^e_d
        end
        ans   = rho_e*m_i_e*m_d_e
    end
    return(ans)
end

#########################
# Rate function for dispersal 
#########################

function rate_d(base,home_ind,away_ind,home_scalar,away_scalar,model)
    ## Error checking
    @assert base > 0 "base rate must be positive" # base rate must be positive 
    @assert home_ind <= num_regions && home_ind > 0 "home region index must be in {1,2,...,num_regions}" # home_ind must be <= num_regions
    @assert away_ind <= num_regions && away_ind > 0 "away region index must be in {1,2,...,num_regions}" # away_ind must be <= num_regions
    @assert home_ind != away_ind "home region index cannot be equal to away region index"
    n_i   = num_species_per_region[home_ind]
    n_j   = num_species_per_region[away_ind]; @assert n_j >= 0 "number of species in away_ind must be >= 0"
    #
    if n_i == 0
        ans = 0 # rate equals 0 if no species in home region
    else
        ## Define rate parameters 
        rho_d      = base # base rate 
        d_d_home   = home_scalar # rate scalar in home region
        d_d_away   = away_scalar # rate scalar in away region
        m_i_d      = rand(Uniform(0.01,1)) # draw density-independent rate randomly [0.01,1]
        ## Compute the overall rates
        if model == "log"
            m_d_d  = (log(n_i)+1)^d_d_home + (log(n_j+1)+1)^d_d_away
        elseif model == "standard"
            m_d_d  = (n_i^d_d_home)*(n_j+1)^d_d_away
        end
        ans   = rho_d*m_i_d*m_d_d
    end
    return(ans)
end

#########################
# Rate function for between-region speciation
#########################

# P.S. the range score function already has a factor of 1/2 to take into account
# left or right children orientation. 

# format for left_range/right_range follows the keys field in dict_keys_b (e.g. [1,0,1] to represent range AC)

function rate_b(base,left_range,right_range,rate_scalar,model)
    ## Error checking
    @assert base > 0 "base rate must be positive" # base rate must be positive 
    @assert left_range + right_range in dict_keys_b "Not a valid split" # left and right ranges must be derived from an ancestral range in dict_keys_b
    @assert left_range != right_range "left and right ranges cannot be equal" # left and right child ranges must be different
    #
    idx_leftrange  = findall(left_range .==1)  # find index of every region that are occupied by the species that forms left range 
    idx_rightrange = findall(right_range .==1) # find index of every region that are occupied by the species that forms right range 
    # @assert all(x->x>0,num_species_per_region[idx_leftrange]) "no species with that given left range" # the number of species with that left range must be > 0
    # @assert all(x->x>0,num_species_per_region[idx_rightrange]) "no species with that given right range" # the number of species with that right range must be > 0
    #
    if any(num_species_per_region[idx_leftrange].==0) == true || any(num_species_per_region[idx_rightrange].==0) == true 
        ans = 0 # if no species in the regions that form either left range or right range, rate = 0
    else
        ## Define rate parameters 
        rho_b      = base # base rate 
        b_d        = rate_scalar # rate scalar 
        m_i_b      = rand(Uniform(0.01,1)) # draw density-independent rate randomly [0.01,1]
        # m_i_b      = 1 # draw density-independent rate randomly [0.01,1]
        ## Compute the range split score f
        f = 0 # initialize the range score for left_range|right_range split 
        for i in 1:length(idx_leftrange)
            for j in 1:length(idx_rightrange)
                if model == "log"
                    m_b_d_ij = (log(num_species_per_region[idx_leftrange[i]]*num_species_per_region[idx_rightrange[j]])+1)^b_d
                    m_b_d_ji = (log(num_species_per_region[idx_leftrange[j]]*num_species_per_region[idx_rightrange[i]])+1)^b_d
                elseif model == "standard"
                    m_b_d_ij = (num_species_per_region[idx_leftrange[i]]+num_species_per_region[idx_rightrange[j]]-1)^b_d
                    m_b_d_ji = (num_species_per_region[idx_leftrange[j]]+num_species_per_region[idx_rightrange[i]]-1)^b_d
                end
                f += 2/(m_b_d_ij + m_b_d_ji) # update f
            end
        end
        f = 1/f
        ## Compute the overall rates
        ans   = rho_b*m_i_b*f
    end
    return(ans)
end

# # Example of non-valid split
# left_range  = [1,1,0,0]
# right_range = [0,1,0,1]
# # Example of valid split
# left_range  = [1,1,0,0]
# right_range = [0,0,1,1]
# num_species_per_region = [2,3,4,5]
# test = rate_b(1,left_range,right_range,1,"log") # should be equal with rate_b(1,right_range,left_range,1)
# println(test)
# xxxx

###########################################
# # Simulate events on tree
###########################################

do_simulation                      = true; # Initial flag to continue simulating events 
start_time                         = 0.0; # Start time of the simulation
start_range                        = floor(Int64,rand(Uniform(1,num_ranges))) # Starting range state at root node 
num_alive_lineage                  = 1; # Start with one species (stem branch)

# base rates 
rho_w = rand(Uniform(0.01,1));
rho_e = rand(Uniform(0.01,1));
rho_d = rand(Uniform(0.01,1));
rho_b = rand(Uniform(0.01,1));

# rate scalars
w_d      = rand(Uniform(-1,1))
e_d      = rand(Uniform(-1,1))
d_d_home = rand(Uniform(-1,1))
d_d_away = rand(Uniform(-1,1))
b_d      = rand(Uniform(-1,1))

# model
model    = "log"

# At root node
num_species_per_range[start_range]      = 1;
range_key_root                          = [k for (k,v) in ranges_inv if v==start_range]; #range key that correspond to range label at the root
reg_idx_root                            = findall(range_key_root[1] .==1 ); # find the region indexes that correspond to that range at the root 
num_species_per_region[reg_idx_root]    .= 1; # update the number of species currently alive in each region 
# println(start_range)
# println(range_key_root)
# println(reg_idx_root)
# println(num_species_per_region)

## Aug 20: Tested loop for computing total sum of rates across all ranges. 
curr_time = start_time; # Initialize current time of simulation 
# while do_simulation == true 
    #############################################
    # # get total sum of rates given the number of currently alive lineages
    #############################################
    w_total = 0 # initialize total sum of rates from within-region speciation
    e_total = 0 # initialize total sum of rates from extinction
    d_total = 0 # initialize total sum of rates from dispersal
    b_total = 0 # initialize total sum of rates from between-region speciation
    for i in 1:num_ranges # loop over each range (total sum of rates = N_A*r_A + ...+ N_ABC*r_ABC)
        range_key = [k for (k,v) in ranges_inv if v==i] # find range that correspond to range label i e.g. [[1,0,0]] for i = 1
        reg_idx   = findall(range_key[1] .==1 ) # find region indexes correspond to that range e.g. reg_idx = 1 for range [1,0,0]
        #
        reg_away  = setdiff([1:1:num_regions;],reg_idx) # find all possible regions for dispersing for species with range i
        for j in 1:length(reg_idx) # for region-specific rates
            global w_total += rate_w(rho_w,reg_idx[j],w_d,model) # e.g. for range {A,B}, w_{A,B} = w_A + w_B
            global e_total += rate_e(rho_e,reg_idx[j],e_d,model) # e.g. for range {A,B}, e_{A,B} = e_A + e_B 
            #
            if reg_away != Int64[] # if the species in range i is not already widespread in all regions 
                for k in 1:length(reg_away)
                    global d_total += rate_d(rho_d,reg_idx[j],reg_away[k],d_d_home,d_d_away,model) #e.g. for range {A,B}, d_AB_j = d_A_C + d_B_C
                end
            else
                global d_total += 0  # if already widepread in all regions
            end
        end
        #
        if count(==(1),range_key[1])==1 # e.g. if range_key = [[1,0,0]], then can't have between-region speciation
            global b_total += 0
        else
            num_left_splits = length(evt_b[range_key[1]])
            for k in 1:num_left_splits # e.g. for range {A,B,C}, can have left range = {A,B} or left_range = {A}
                left_range  = evt_b[range_key[1]][k] 
                right_range = range_key[1] - left_range
                global b_total += 2*rate_b(rho_b,left_range,right_range,b_d,model) # 2 is there because {A}|{B} == {B}|{A}
            end
        end
    end
    total_rate = sum(num_species_per_range)*(w_total+e_total+d_total+b_total) # total sum of rates across ranges 
    # println(w_total)
    # println(e_total)
    # println(d_total)
    # println(b_total)
    # xxxx

    #############################################
    # # draw time to next event
    #############################################
    wait_time = rand(Exponential(1/total_rate),1)[1] # draw waiting time till a next event 
    global curr_time += wait_time # update current time 
    # println(wait_time)

    # # draw next event type
    prob_w   = w_total/total_rate # probability that the next event is within-region speciation
    prob_e   = e_total/total_rate # probability that the next event is extinction 
    prob_d   = d_total/total_rate # probability that the next event is dispersal 
    prob_b   = b_total/total_rate # probability that the next event is between-region speciation 
    #
    rand_num = 10 # some random big number to initiate the while loop 
    while all([prob_w,prob_e,prob_d,prob_b].<rand_num) == true # keep re-drawing the random number until we get a valid one. Is this correct?
        global rand_num = rand(Uniform(0,1)) # draw random number within [0,1] uniformly 
    end
    #
    if rand_num <= prob_w && prob_w == max(prob_w,prob_e,prob_d,prob_b) && any(num_species_per_range.!=0) == true # if next event is within-reg spec and there's alive species with some ranges
        curr_event = "within-region" # assign next event to be within-region speciation 
    elseif rand_num <= prob_e  && prob_e == max(prob_w,prob_e,prob_d,prob_b) && any(num_species_per_range.!=0) == true # if next event is extinction and there's alive species with some ranges
        curr_event = "extinction" # assign next event to be extinction
    elseif rand_num <= prob_d  && prob_d == max(prob_w,prob_e,prob_d,prob_b) && any(num_species_per_range[1:end-1].!=0) == true # check not including the last entry in num_species_per_range because that species cannot disperse
        curr_event = "dispersal" # assign the next event to be dispersal
    elseif rand_num <= prob_b  && prob_b == max(prob_w,prob_e,prob_d,prob_b) && any(num_species_per_range[num_regions+1:end].!=0) == true # check not including the first num_regions entries because species have range size 1 (can't speciate between region)
        curr_event = "between-region" # assign the next event to be between-region 
    end
    # # Aug 21: works as intended 
    # println("random number = ",rand_num)
    # println("prob w = ",prob_w)
    # println("prob e = ",prob_e)
    # println("prob d = ",prob_d)
    # println("prob b = ",prob_b)
    # println("num species range = ",num_species_per_range)
    # println("current event = ",curr_event)

    #############################################
    # # assign the next event type to a particular alive lineage
    #############################################

    # if curr_event == "within-region"
    #     # ..... choose which lineage
    # elseif curr_event == "extinction"
    #     # .... choose which lineage
    #     if size == 1
    #         curr_event == "lineage extinction"
    #     elseif size != 1
    #         curr_event == "local extinction"
    #     end
    # elseif curr_event == "dispersal"
    #     # .... choose which lineage
    # elseif curr_event == "between-region"
    #     # .... choose which lineage
    # end
    
#     # update current simulation state and record parent-child-event_type (use dataframe?)
    
#     # check if we stop simulation
#     if stop_condition == "time" and curr_time > max_time  # if condition on max_time
#         do_simulation = false
#     elseif stop_condition == "taxa" and curr_taxa > max_taxa # if condition on max_taxa
#         do_simulation = false
#     end
    
#     # stop simulation if condition met
#     if do_simulation == false 
#         break
#     end
    
# end

# # output

# # write tree to file

# # write ranges to file

# # write parameters to file

# # write summary statistics to file