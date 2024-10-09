# This script will simulate phylogenetic datasets under a 
# density-dependent GeoSSE model.

# This Julia script already has while loop (so we need remove some global assignment since it's not inside a loop
# e.g. w_total, e_total etc)

# still has a bug where the next event is picked from a state that has no species in it (from num_species_per_range)
# but its  state_event_tot_rate  has non-zero for that corresponding entry involving that range
# (the issue is from lookup table for d when there is already a recorded entry, so I turn it off for now and it fixed)

# and fixed a bug in computing rate_b (in m_b_d_ij  and m_b_d_ji out of range)

# this script sim_test_v5.jl 
# - added @time macro to measure script's run time 
# - added condition if curr_time > max_time then we do not draw the next event
# - increased size of lookup_w and lookup_e to accommodate more species in a region
# - added condition on max_taxa 

# this script sim_test_v6.jl
# - added 3 more columns to history_df (branch_lengths, child_1_idx and child_2_idx)
# - relabeled column "state" to "final_state" for more clarity

# this script sim_test_v7.jl
# - convert history_df to newick string and save the string to a file system 

# this script sim_test_v8.jl
# - added command line argument to generate tree with index i

# this script sim_test_v9.jl
# - added printing line to indicate that tree i has been simulated 
# - added command line argument 
# - added option to re-simulate (e.g. in case of extinction)
# - added printing line to indicate if tree i is being re-simulated 
# - write tip state data to csv file
# - write other summary statistics and rate parameter to files 
# - added printing line for base rate, rate scalars, current time, and num_alive_lineage for debugging 
# Note: for some choice of rate parameters, events occcur in such short amount of time (trees grows very large)

# this script sim_test_v10.jl
# -. only draw negative rate scalars to prevent runaway in taxa count.
# -. added continuous time points from timebin and max_time to record number of species across regions in cts time 
# -. added continuous time points from timebin and max_time to record number of species across ranges in cts time 
# -. added event times to region_count_history and range_count_history
# -. lower base extinction rate to avoid trees keep becoming extinct under a certain set of base event rates
# -. fixed typo on m_d_d (changed + to *)
# -. added the correct extensions for the output i.e. labels.csv, dat.csv and .csv 
# -. reduced the max runtime from 50 to 35

# this script sim_test_v11.jl
# -. increased n in lookup_w and lookup_e from 10000 to 100000
# -. added new column to history_df for event_type at each cladogenetic event 
# -. write cladogenetic event history over event times (from history_df)
# -. fixed data frame structures before writing to csv (i.e. remove any, and have separate columns)
# -. allow event rates to accept base rate of 0 (not just positive base rates)
# -. extend the runtime in order to reach the local equilibrium capacity under the testing described below 
# -. fixed root state to always start in range A for testing purpose described below:
# In this, we valide n_i(t) for i = A under scenario of no dispersal and between-region, and constant extinction rate
# we find the base rate for w event that corresponds to n_i = 100 in the limit. 

# TBA: maximum number of times tree i can be resimulated? (i.e. so we do not have an infinite loop in case that set of parameters is impossible to generate a valid tree)

# @time begin
    #####################
    # Load julia packages
    #####################
    using Combinatorics # Package for combinatorics/permutations
    using OrderedCollections # Package for ordering dictionary 
    using DataStructures # Package for preserving dictionary order by its insertion order 
    using Distributions # Package for using probability distributions 
    using DataFrames # Package for dataframe
    using CSV # Package to write to csv
    
    ####################
    # Command line arguments 
    ####################
    args = ARGS
    sim_id     = map(x->string(x), args)[1] # tree id
    resimulate = args[2] # whether to resimulate in case of tree extinction or not
    timebin    = args[3] # time bin for tracking number of species across region/range in continuous time 
    # sim_id = string(1)
    # resimulate = true
    # timebin = 10

    #############
    # File system
    #############
    # base_dir = "~/projects/density_geosse";
    base_dir = "/Users/albertsoewongsono/Documents/Code\\ Testing/density_geosse/code/jl"
    out_dir  = base_dir*"output/";
    ################
    # Initial setting
    ################
    min_taxa       = 5; # Minimum number of extant taxa to simulate 
    # max_taxa       = 100; # Maximum number of extant taxa to simulate
    max_taxa       = 10; # Maximum number of extant taxa to simulate
    # max_time       = 0.5; # Maximum simulation time 
    # max_time       = 35; # Maximum simulation time 
    max_time       = 500; # Maximum simulation time 
    stop_condition = "time" # If time is the stopping condition
    # stop_condition = "numtaxa" # If the final number of extant taxa is the stopping condition
    
    #############
    # Model setup
    #############
    # Define state space 
    # num_regions = 4; # Number of regions 
    num_regions = 3; # Number of regions 
    num_ranges  = 2^num_regions - 1; # Number of possible ranges given num_regions
    num_species_per_range  = fill(0, num_ranges); # Empty vector to store number of species in each range
    # num_species_per_range = [10,0,0,0,0,0,0]
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
    #     (0,0,1):3,
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
    # println(evt_w)
    # xxxx
    
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
        else # if the ancestral species is widespread, can speciate
            anc_range_size = count(==(1),dict_keys_b[k]) # ancestral range size for index k in key field
            # looping for each different size of child range after speciation.
            # Divide by 2 because the range size of the other child can be derived from the range of the other one. 
            idx_anc = findall(dict_keys_b[k].== 1 ); # find locations of 1's in ancestral range k
            #
            child_range_b = Vector{Int}[]; # initialize vector of vectors to store all child ranges from ancestral range k 
            for l in 1:floor.(Int64,anc_range_size/2) # loop over each possible range split size for the child range e.g. for k = 4, we can have split of size 1&2 (not including its compelement of size 3)
                    child_init = zeros(Int64,num_regions); # Initialize vector of zeros 
                    child_init[1:l] .= 1; # Assign 1's equal to the range split size 
                    child_comb = unique(permutations(child_init)); # permute all possibilities of combinations of 1's for that split size 
                    if l != anc_range_size - l # case: not equal split between left and right child range size e.g. l|m where l!=m
                        for m in 1:length(child_comb) # loop over this permutation 
                            idx_child = findall(child_comb[m].== 1 );
                            test_idx = [x in idx_anc for x = idx_child] # Check if that split of size l is valid (i.e. if its contained in ancestral range)
                            if all(test_idx) == true 
                                push!(child_range_b,child_comb[m]) # If yes, then we accept that split of size l as a valid split
                            end
                        end
                    else # case: equal split between left and right child range size e.g. l|m where l = m 
                        track_range = Vector{Int}[] # initialize empty vector of vectors 
                            for n in 1:length(child_comb) # loop over this permutation 
                                idx_child = findall(child_comb[n].== 1 );
                                test_idx = [x in idx_anc for x = idx_child] # Check if that split of size l is valid (i.e. if its contained in ancestral range)
                                if all(test_idx) == true 
                                    push!(child_range_b,child_comb[n]) # If yes, then we accept that split of size l as a valid split
                                    push!(track_range,child_comb[n]) # If yes, add that to track_range
                                    break
                                end
                            end
                        child_comb_reduce = setdiff(child_comb,track_range) # remove the newly added child range from the pool in child_comb
                        for p in 1:length(child_comb_reduce) # loop over the remaining
                            idx_child = findall(child_comb_reduce[p].== 1 );
                            test_idx = [x in idx_anc for x = idx_child] # Check if that split of size l is valid (i.e. if its contained in ancestral range)
                            if (all(test_idx) && child_comb_reduce[p] ∉ [dict_keys_b[k]].-track_range) == true
                                push!(child_range_b,child_comb_reduce[p])
                                push!(track_range,child_comb_reduce[p])
                            end
                        end
                    end
            end
            evt_b[dict_keys_b[k]] = child_range_b; # append keys and value fields to dictionary 
        end
    end 
    
    # println(evt_b)# Note for range size 4, the complement (range size 2) is still recorded 
    # print(dict_keys_b)
    # println(evt_b[[1,1,1,0]])
    # println(typeof(evt_b[[1,1,1,0]]))
    # println([1,1,1,0]-evt_b[[1,1,1,0]][1])
    # println(length(evt_b)) 
    # xxx
    
    ###########################################
    # # Define rates for event space in each model variant 
    ###########################################
    
    # log-normal distribution (mean = m, std = std) in Julia 
    # p.s. use lognormal to ensure rate is always positive
    function myLogNormal(m,std)
        γ = 1+std^2/m^2
        μ = log(m/sqrt(γ))
        σ = sqrt(log(γ))
        return LogNormal(μ,σ)
    end
    
    # println(rand(myLogNormal(1,1)))
    
    #########################
    # Rate function for within-region speciation 
    #########################
    
    # note: use lognormal distribution with mean 1 (for no effect from non-density)
    # use flag (true or false) for using density-independent or not (rate == 1)
    
    function rate_w(base,region_ind,species_count_ind,rate_scalar,only_density,model)
        ## Error checking
        # @assert base > 0 "base rate must be positive" # base rate must be positive 
        @assert base >= 0 "base rate must be positive" # base rate must be positive 
        @assert region_ind <= num_regions && region_ind > 0 "region index must be in {1,2,...,num_regions}" # region ind must be <= num_regions
        n_i   = species_count_ind # Number of species currently in region_ind
        if n_i == 0
            ans = 0 # rate equals 0 if no species 
        else
            ## Define rate parameters 
            rho_w = base # base rate 
            w_d   = rate_scalar # rate scalar for within-region speciation
            if only_density == true 
                m_i_w = 1 # means no density-independent factors 
            elseif only_density == false
                m_i_w = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 
            end 
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
    # test = rate_w(0.1,1,0,0.1,false,"log") # 0 species in region 1
    # test2 = rate_w(0.1,1,10,0.1,true,"log") # 10 species in region 1
    # println(test)
    # println(test2)
    # xxx
    
    #########################
    # Rate function for extinction
    #########################
    
    function rate_e(base,region_ind,species_count_ind,rate_scalar,only_density,model)
        ## Error checking
        # @assert base > 0 "base rate must be positive" # base rate must be positive
        @assert base >= 0 "base rate must be positive" # base rate must be positive 
        @assert region_ind <= num_regions && region_ind > 0 "region index must be in {1,2,...,num_regions}" # region ind must be <= num_regions
        n_i   = species_count_ind # Number of species currently in region_ind
        if n_i == 0
            ans = 0 # rate equals 0 if no species 
        else # if there is species 
            ## Define rate parameters 
            rho_e = base # base rate 
            e_d   = rate_scalar # rate scalar for extinction
            if only_density == true 
                m_i_e = 1 # means no density-independent factors 
            elseif only_density == false
                m_i_e = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 
            end 
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
    
    function rate_d(base,home_ind,away_ind,home_scalar,away_scalar,only_density,model)
        ## Error checking
        # @assert base > 0 "base rate must be positive" # base rate must be positive 
        @assert base >= 0 "base rate must be positive" # base rate must be positive 
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
            if only_density == true 
                m_i_d = 1 # means no density-independent factors 
            elseif only_density == false
                m_i_d = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 
            end 
            ## Compute the overall rates
            if model == "log"
                m_d_d  = (log(n_i)+1)^d_d_home * (log(n_j+1)+1)^d_d_away
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
    function rate_b(base,left_range,right_range,rate_scalar,only_density,model)
        ## Error checking
        # @assert base > 0 "base rate must be positive" # base rate must be positive 
        @assert base >= 0 "base rate must be positive" # base rate must be positive 
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
            if only_density == true 
                m_i_b = 1 # means no density-independent factors 
            elseif only_density == false
                m_i_b = rand(myLogNormal(1,1)) # draw density-independent rate from lognormal with mean 1, std 1 
            end 
            # m_i_b      = 1 # draw density-independent rate randomly [0.01,1]
            ## Compute the range split score f
            f = 0 # initialize the range score for left_range|right_range split 
            for i in 1:length(idx_leftrange)
                for j in 1:length(idx_rightrange)
                    if model == "log"
                        m_b_d_ij = (log(num_species_per_region[idx_leftrange[i]]*num_species_per_region[idx_rightrange[j]])+1)^b_d
                        m_b_d_ji = (log(num_species_per_region[idx_rightrange[j]]*num_species_per_region[idx_leftrange[i]])+1)^b_d
                    elseif model == "standard"
                        m_b_d_ij = (num_species_per_region[idx_leftrange[i]]+num_species_per_region[idx_rightrange[j]]-1)^b_d
                        m_b_d_ji = (num_species_per_region[idx_rightrange[j]]+num_species_per_region[idx_leftrange[i]]-1)^b_d
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
    # println(rate_b(1,[1,0,0],[0,1,1],1,true,"log"))
    # xxx
    # # Example of non-valid split
    # left_range  = [1,1,0,0]
    # right_range = [0,1,0,1]
    # # Example of valid split
    # left_range  = [1,1,0,0]
    # right_range = [0,0,1,1]
    # num_species_per_region = [2,3,4,5]
    # test = rate_b(1,left_range,right_range,1,true,"log") # should be equal with rate_b(1,right_range,left_range,1)
    # println(test)
    # xxxx
    
    ###########################################
    # # Simulate events on tree (set up)
    ###########################################
    
    do_simulation                      = true; # Initial flag to continue simulating events 
    start_time                         = 0.0; # Start time of the simulation (0 = past time)
    # start_range                        = floor(Int64,rand(Uniform(1,num_ranges))); # Starting range state at root node e.g. 1,2,...,7 for 3 region
    start_range                        = 1;
    num_alive_lineage                  = 1; # Start with one species (stem branch)
    
    ### sim1
    # base rates 
    # rho_w = rand(Uniform(0.01,1));
    # rho_e = rand(Uniform(0.01,0.1));
    # rho_d = rand(Uniform(0.01,1));
    # rho_b = rand(Uniform(0.01,1));

    ### sim2 (no dispersal, no between-region speciation, constant extinction rate)
    rho_w = 3.141793;
    rho_e = 0.1;
    rho_d = 0;
    rho_b = 0;

    # println("base w ",rho_w)
    # println("base e ",rho_e)
    # println("base d ",rho_d)
    # println("base b ",rho_b)

    ### sim1
    # rate scalars (draw from standard normal distribution)
    # only consider negative direction to avoid runaway in taxa count
    # p.s. rate scalar == 0 -> no density-dependent effect 
    # w_d      = rand(Normal(-4,1));
    # e_d      = rand(Normal(-4,1));
    # d_d_home = rand(Normal(-4,1));
    # d_d_away = rand(Normal(-4,1));
    # b_d      = rand(Normal(-4,1));

    ### sim2 (no dispersal, no between-region speciation, constant extinction rate)
    w_d      = -2;
    e_d      = 0;
    d_d_home = rand(Normal(-4,1));
    d_d_away = rand(Normal(-4,1));
    b_d      = rand(Normal(-4,1));

    # println("scalar w ",w_d)
    # println("scalar e ",e_d)
    # println("scalar d_i ",d_d_home)
    # println("scalar d_j ",d_d_away)
    # println("scalar b ",b_d)
    
    # model
    model        = "log";
    # only_density = false; # not including other density-independent factors when computing rates 
    only_density = true; # just for testing the output, for easier computation since density_indep rates = 1
    
    # set up lookup table for  pre-computed region-based within-region speciation rates (only need to set up for one region)
    lookup_w = Float64[]; # initialize empty vector to store rate_w
    for i in 1:100001 # loop over each possible number of species in a region (from n_i = 0 to n_i = 100000)
        value_w = rate_w(rho_w,1,i-1,w_d,only_density,model);
        push!(lookup_w,value_w); # append to lookup_w
    end
    # println(lookup_w[1:5])
    # xxx
    
    # set up lookup table for  pre-computed region-based extinction rates (only need to set up for one region)
    lookup_e = Float64[]; # initialize empty vector to store rate_e
    for i in 1:100001 # loop over each possible number of species in a region (from n_i = 0 to n_i = 100000)
        value_e = rate_e(rho_e,1,i-1,e_d,only_density,model);
        push!(lookup_e,value_e); # append to lookup_e
    end
    # println(lookup_e[1:5])
    # xxx
    
    # initialize lookup table for dispersal rates (cached during simulation and revisit the combination if ever)
    # e.g. combination in num_species_per_region
    # only need to set up for number of species in any home location i and number of species in any away location j
    # New: this is already multiplied by number of species in each range e.g. if i = {AB}, then the two rows corresponds
    # to N_AB*d_AC and N_AB*d_BC
    
    lookup_d = DataFrame(range_label = Int64[],num_species_i = Int64[],num_species_j = Int64[],relative_dispersal_ij = Float64[]);
    
    # initialize lookup table for between-region speciation rates (cached during simulation and revisit the combination if ever)
    # p.s. (left_range,right_range) and (right_range,left_range) will have the same rate
    nodenames = ["ancestral_range","left_range","right_range"]
    for i in 1:num_regions
        reg_name = "num_spec_region_$i";
        push!(nodenames,reg_name)
    end
    push!(nodenames,"relative between-reg rate")
    lookup_b = DataFrame([name => [] for name in nodenames])
    #
    for j in (num_regions+1):num_ranges #loop over ancestral range with size > 1
        for k in 1:length(evt_b[dict_keys_b[j]]) # loop over each left child species' range for the given ancestral range
            new_row = ([dict_keys_b[j],evt_b[dict_keys_b[j]][k],dict_keys_b[j]- evt_b[dict_keys_b[j]][k]]) # for a given left child range, find its right species' range 
            for l in 1: num_regions + 1
                new_row = vcat(new_row,0) # assign placeholder for number of species in each region and rate to be equal to 0 (temporary)
            end
            push!(lookup_b,new_row) # push to the look up table
        end
    end
    # println(lookup_b)
    # xx
    
    ##########################
    # Update event at root node
    ##########################
    num_species_per_range[start_range]      = 1;
    # num_species_per_range = [1,0,0,1,0,0,0]; # for testing total rate is computed correctly
    range_key_root                          = [k for (k,v) in ranges_inv if v==start_range]; #range key that correspond to range label at the root e.g. (1,1,0) for label 4
    reg_idx_root                            = findall(range_key_root[1] .==1 ); # find the region indexes that correspond to that range at the root e.g. region_ind = 1,2 for range (1,1,0)
    num_species_per_region[reg_idx_root]    .= 1; # update the number of species currently alive in each region e.g. 1 species with range {A,B} means 1 species in A and 1 species in B
    # num_species_per_region = [2,1,0]; # for testing total rate is computed correctly
    # println(start_range)
    # println(range_key_root)
    # println(reg_idx_root)
    # println(num_species_per_region)
    
    #
    range_count_history  = [] # record num_species_per_range across time  
    region_count_history = [] # record num_species_per_region across time 
    #

    push!(range_count_history,deepcopy(num_species_per_range))
    range_count_history[1] = insert!(Float64[x for x in range_count_history[1]],1,start_time)

    push!(region_count_history,deepcopy(num_species_per_region))
    region_count_history[1] = insert!(Float64[x for x in region_count_history[1]],1,start_time)

    # create array of vectors containing information about num_species across regions in cts time
    timebin                 = deepcopy(parse(Int,timebin)) # convert string (from args) to int
    cts_time_pts            = collect(range(0,max_time,length=timebin)) # continuos time points according to timebin

    species_cts_time_region = [] # initialize array of vectors containing information about num_species across regions in cts time
    for i in 1:timebin
        added_vec     = fill(missing,num_regions)
        added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
        insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
        #
        push!(species_cts_time_region,added_vec_mod)
    end 
    species_cts_time_region[1][2:(num_regions+1)]=deepcopy(num_species_per_region)
    #
    species_cts_time_range = [] # initialize array of vectors containing information about num_species across ranges in cts time
    for i in 1:timebin
        added_vec     = fill(missing,num_ranges)
        added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
        insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
        #
        push!(species_cts_time_range,added_vec_mod)
    end 
    species_cts_time_range[1][2:(num_ranges+1)]=deepcopy(num_species_per_range)
   
    ## Aug 20: Tested loop for computing total sum of rates across all ranges. 
    ## Aug 23: Tested loop for computing w_total,e_total,d_total working, using the lookup table
    
    curr_time = start_time; # Initialize current time of simulation p.s. time starts from 0 (past)
    ###############################################
    # Dataframe to record parent-child relationships 
    ###############################################
    # idx = branch label, time_start = starting time of the branch, time_end = end time of the branch (missing if still continues)
    # final_state = range state on the branch just before an event on that branch, extinct_status = whether its already extinct (1) or not (0) and "missing" if no event on the branch yet
    # split_state = whether it corresponds to a splitting event (1) or not (0) and "missing" if no event occuring on the branch yet,
    # branch_length = length of branch label idx (time_end - time_start), child_1_idx = idx of 1st child (missing if no speciation), child_2_idx = idx of 2nd child (missing if no speciation)
    history_df = DataFrame(idx = [],time_start = [],time_end = [],final_state = [],extinct_status = [],split_status = [],branch_length = [],child_1_idx = [],child_2_idx = [],event_type = []) 
    push!(history_df,(1,start_time,missing,start_range,missing,missing,missing,missing,missing,missing)) # add first event at the root node (root node)
    # history_df[1,:time_end]=1   
    # println("example df = ",history_df)
    # xxx
    
    ###############################################
    # Dataframe to record anagenetic transitions (local extinction on widespread range or dispersal) 
    ###############################################
    # idx = branch label associated with the event, time = time when the event occurs 
    # from_state = range state before the event, to_state = range state after the event 
    anagen_df = DataFrame(idx = [],time = [], from_state = [], to_state = [])
    
    while do_simulation == true 
        #############################################
        # # get total sum of rates given the number of currently alive lineages
        #############################################
        w_total = 0 # initialize total sum of rates from within-region speciation
        e_total = 0 # initialize total sum of rates from extinction
        d_total = 0 # initialize total sum of rates from dispersal
        b_total = 0 # initialize total sum of rates from between-region speciation
        ############################################
        # # Create vector for state-event relationship 
        ############################################
        # P.S. 1st idx = range idx (e.g. A,B,AB,etc), 2nd idx = event type (w,e,d,b), 3rd idx = left range (for b) or home location for (for d) or region index (which region it occurs for w and e)
        # 4th idx = right range (for b) or away location (for d) or 0 (for w and e)
        state_event_col_name      = Tuple{Any,Any,Any,Any}[] 
        state_event_tot_rate      = zeros(0) # initialize vector to store values from each index in state_event_col_name
        ##
        for i in 1:num_ranges # loop over each range (total sum of rates = N_A*r_A + ...+ N_ABC*r_ABC)
            range_key = [k for (k,v) in ranges_inv if v==i] # find range that correspond to range label i e.g. [[1,0,0]] for i = 1
            reg_idx   = findall(range_key[1] .==1 ) # find region indexes correspond to that range e.g. reg_idx = 1 for range [1,0,0]
            #
            reg_away  = setdiff([1:1:num_regions;],reg_idx) # find all possible regions for dispersing for species with range i
            for j in 1:length(reg_idx) # for region-specific rates
                n_j = num_species_per_region[reg_idx[j]] # find the number of currently living species in region reg_idx[j]
                # compute N_i * w_j and N_i * e_j where i = range state and j in i  
                prod_ni_wj = num_species_per_range[i]*lookup_w[n_j+1] # e.g. for range {A,B}, w_{A,B} = n_AB * (w_A + w_B) with n_j species in A or B from lookup table
                prod_ni_ej = num_species_per_range[i]*lookup_e[n_j+1] # e.g. for range {A,B}, e_{A,B} = n_AB * (e_A + e_B) with n_j species in A or B from lookup table
                w_total += prod_ni_wj # e.g. w_total = N_A*(w_A) + N_B*(w_B) + N_AB*(w_A+w_B)
                e_total += prod_ni_ej # e.g. e_total = N_A*(e_A) + N_B*(e_B) + N_AB*(e_A+e_B)
                # append each state-event type to state_event, along with its column name
                push!(state_event_col_name,(i,"w",reg_idx[j],0),(i,"e",reg_idx[j],0)) # e.g. (1,"w",1,0) denotes N_A * w_A or (2,"e",2,0) denotes N_B * w_B
                append!(state_event_tot_rate,prod_ni_wj,prod_ni_ej)
                #
                if reg_away != Int64[] # if the species in range i is not already widespread in all regions 
                    for k in 1:length(reg_away)
                        # check if the given combination of n_j and n_k already has its pre-computed dispersal rate in lookup_d
                        idx_match_d = (lookup_d.range_label .== i) .& (lookup_d.num_species_i .== num_species_per_region[reg_idx[j]]) .& (lookup_d.num_species_j .== num_species_per_region[reg_away[k]])
                        # if  all(!,idx_match_d) # if that n_j and n_k does not exist yet, compute the rate and store in lookup_d
                            disp_rate = rate_d(rho_d,reg_idx[j],reg_away[k],d_d_home,d_d_away,only_density,model) # compute dispersal rate from home to away
                            # compute N_i * d_jk where i = range state, and j in i, k not in i
                            prod_ni_dik = num_species_per_range[i]*disp_rate 
                            d_total += prod_ni_dik  # e.g. for range {A,B}, d_AB_j = n_AB*(d_A_C + d_B_C) or d_A_j = n_A*(d_A_B + d_A_C)
                            push!(lookup_d,[i,num_species_per_region[reg_idx[j]],num_species_per_region[reg_away[k]],disp_rate])
                            # append each state-event type to state_event, along with its column name
                            push!(state_event_col_name,(i,"d",reg_idx[j],reg_away[k])) # e.g. (1,"d",1,2) denotes N_A * d_AB or (4,"d",1,3) denotes N_AB * d_AC
                            append!(state_event_tot_rate,prod_ni_dik)
                        # else # if already exist, recall rate from lookup_d
                        #     prod_ni_dik = lookup_d.relative_dispersal_ij[findall(idx_match_d .== 1)][1]
                        #     d_total += prod_ni_dik
                        #     # append each state-event type to state_event, along with its column name
                        #     push!(state_event_col_name,(i,"d",reg_idx[j],reg_away[k]))
                        #     append!(state_event_tot_rate,prod_ni_dik)
                        # end
                    end
                else # if already widepread in all regions
                    prod_ni_dik = num_species_per_range[i]*0  
                    d_total += prod_ni_dik
                    # Then we do not append to state_event for that range state and that particular event type 
                end
            end
            #
            if count(==(1),range_key[1])==1 # e.g. if range_key = [[1,0,0]], then can't have between-region speciation
                # compute N_i * b_leftrange_right_range where leftrange U rightrange = i
                prod_ni_b =  num_species_per_range[i]*0
                b_total += prod_ni_b
                # Then we do not append to state_event 
            else
                num_left_splits = length(evt_b[range_key[1]])
                for k in 1:num_left_splits # e.g. for range {A,B,C}, can have left range = {A,B} or left_range = {A} (but not considering their change in orientation)
                    left_range  = evt_b[range_key[1]][k] 
                    right_range = range_key[1] - left_range
                    # # Note: 26 Aug, have not added option to push to lookup table for future references
                    # # idx_match_b = (lookup_b.ancestral_range .== range_key[1]) .& (lookup_b.left_range .== left_range) .& (lookup_b.right_range .== right_range)
                    #
                    # compute N_i * b_leftrange_right_range where leftrange U rightrange = i
                    prod_ni_b = num_species_per_range[i]*rate_b(rho_b,left_range,right_range,b_d,only_density,model)
                    b_total +=  prod_ni_b # e.g. i = {A,B}, then its N_AB*(b_A_B).agreed: no left-right orientation considered.
                    # append each state-event type to state_event, along with its column name
                    left_range_id = [v for (k,v) in ranges_inv if k==left_range] # e.g. 1 for range [1,0,0]
                    right_range_id = [v for (k,v) in ranges_inv if k==right_range]
                    push!(state_event_col_name,(i,"b",left_range_id[1],right_range_id[1])) # e.g. (4,"b",1,2) denotes N_AB * b_A_B
                    append!(state_event_tot_rate,prod_ni_b)
                    # # # push to the lookup table
                    # # added_row = [range_key,left_range,right_range]
                    # # for num_r in 1:num_regions
                    #     # added_row = vcat(added_row,num_species_per_region[num_r])
                    # # end
                    # # push!(lookup_b,added_row)
                end
            end
        end
        # total sum of rates across ranges 
        # e.g. w_total = N_A*w_A + N_B*w_B + N_AB*(w_A+w_B) (loop over each range, and all event associated with each range)
        # e.g. e_total = N_A*e_A + N_B*e_B + N_AB*(e_A+e_B)
        # e.g. d_total = N_A*(d_AB) + N_B*(d_BA)
        # e.g. 
        total_rate = w_total+e_total+d_total+b_total # Note: this should match with sum(state_event_tot_rate) - checked 
        # Check these result total_rate and individual component by hand (Aug 28. Checked for correctness)
        # println("w_total = ",w_total)
        # println("e_total = ",e_total)
        # println("d_total = ",d_total)
        # println("b_total = ",b_total)
        # println("range_vec = ",num_species_per_range)
        # println("region_vec = ",num_species_per_region)
        # println("base w = ",rho_w) 
        # println("base e = ",rho_e) 
        # println("base d = ",rho_d) 
        # println("base b = ",rho_b) 
        # println("rate scalar w = ",w_d)      
        # println("rate scalar e = ",e_d)      
        # println("d_d_home = ",d_d_home) 
        # println("d_d_away = ",d_d_away) 
        # println("between = ",b_d)      
        # println("vec_name",state_event_col_name)
        # println("length_vec_name = ",length(state_event_col_name))
        # println("vec_value",state_event_tot_rate)
        # println("length_vec_value = ",length(state_event_tot_rate))
        # println(total_rate)
        # println(sum(state_event_tot_rate))
        # xxxx
    
        #############################################
        # # draw time to next event
        #############################################
        wait_time = rand(Exponential(1/total_rate),1)[1] # draw waiting time till a next event 
        global curr_time += wait_time # update current time 
        # println("current time ", curr_time)
        # println(wait_time)
        #
        #############################################
        # # sampling next event type and associated state at the updated curr_time 
        #############################################
        # PS. each cell in state_event_tot_prob corresponds with a product of N_i and an event rate associated with the range i.
        if (curr_time <= max_time && stop_condition == "time") || (stop_condition == "numtaxa")# If conditioned on time, and curr_time > max_time, do not draw next event. 
            state_event_tot_prob  =  state_event_tot_rate/sum(state_event_tot_rate) # convert rates to probabilities (e.g. N_A*w_A | N_AB*w_A | ....| ...)
            rand_u = rand(Uniform(0,1)) # draw random number between [0,1]
            next_event = Tuple{Any,Any,Any}[] 
            next_state = zeros(Int64,0)
            for move in 1:length(state_event_tot_prob) # sampling from probability distribution
                # global rand_u = rand_u - state_event_tot_prob[move]
                rand_u = rand_u - state_event_tot_prob[move]
                if rand_u <0 # we have found over next event type and state associated with the event type 
                    push!(next_event,state_event_col_name[move][2:4]) # find the next event type
                    append!(next_state,state_event_col_name[move][1]) # find the associated state
                    break
                end
            end 
            # # Aug 30: checked it's working as intended
            # # p.s. for range state and event type that is not possible (e.g. N_A with between-region), this is not recorded
            # # in state_event_col_name and state_event_tot_rate (state_event_tot_rate)
    
            # println("name =", state_event_col_name)
            # println("num species range = ",num_species_per_range)
            # println("next event = ",next_event[1])
            # println("next_state = ",next_state[1])
            # println("event type = ",next_event[1][1])
            # xxx
    
            #############################################
            # # assign the next event type to a particular alive lineage that are available 
            #############################################
            #############################################
            # # Enable these prints if want to check to see lineage history at each iteration 
            ############################################
            # find row indexes from dataframe that satisfy the conditions (no event occuring yet and has the same range state as next_state)
            # println("parent-child history = ",history_df)
            # println("anagenetic history = ", anagen_df)
            # println("num species per range = ",num_species_per_range)
            # println("num species per region",num_species_per_region)
            # println("chosen state = ",next_state[1])
            # println("next event = ",next_event[1])
    
            row_idx    = intersect(findall(ismissing,history_df.extinct_status),findall(ismissing,history_df.split_status),
                                findall(history_df.final_state .== next_state[1])) 
            row_chosen = rand(row_idx) # randomly sampling which branch in history_df[row_idx,:].idx that undergoes the next event 
    
            # # Example
            # push!(history_df,(2,start_time,missing,8,missing,missing))
            # push!(history_df,(3,start_time,missing,8,missing,missing))
            # push!(history_df,(4,start_time,missing,8,0,missing))
            # push!(history_df,(5,start_time,missing,8,missing,0))
            # push!(history_df,(6,start_time,missing,next_state[1],missing,missing))
            # #
            # row_idx = intersect(findall(ismissing,history_df.extinct_status),findall(ismissing,history_df.split_status),findall(history_df.final_state .== next_state[1]))
            # println(history_df)
            # #
            # history_df[row_idx,:extinct_status] .= 1 # replace a column entry from those row indexes with a new value
            # history_df[row_idx,:split_status] .= 0
            # println(history_df)
            # println(next_state[1])
            # println("row index = ",row_idx)
            # println(rand(row_idx))
            # xxx
            # 
    
            if next_event[1][1] == "w" # if the next event is within-region speciation (Sep 5: already checked)
                max_idx = maximum(history_df.idx) # find the maximum idx from current recorded branch labels in history_df
                push!(history_df,(max_idx+1,curr_time,missing,next_state[1],missing,missing,missing,missing,missing,missing)) # child 1 inherits the parent's range
                push!(history_df,(max_idx+2,curr_time,missing,next_event[1][2],missing,missing,missing,missing,missing,missing)) # child 2 inherits the other range in case parent is widespread
                # update the parent branch's status
                # P.S. row_chosen must be transformed into a vector instead of scalar
                history_df[[row_chosen],:time_end]       .= curr_time
                history_df[[row_chosen],:extinct_status] .= 0
                history_df[[row_chosen],:split_status]   .= 1
                history_df[[row_chosen],:branch_length]  .= curr_time - deepcopy(history_df[[row_chosen],:time_start][1])
                history_df[[row_chosen],:child_1_idx]    .= deepcopy(max_idx+1)
                history_df[[row_chosen],:child_2_idx]    .= deepcopy(max_idx+2)
                history_df[[row_chosen],:event_type]     .= "w"
                # println(row_chosen)
                # println(history_df)
                # println(history_df[row_chosen,:time_end])
                # xxx
            elseif next_event[1][1] == "e" && next_state[1] > num_regions # if the next event is extinction on a branch with widespread range (Sep 5: already checked)
                # update the chosen branch's status 
                original_range  = [k for (k,v) in ranges_inv if v==next_state[1]] # the original range of the branch before the extinction e.g. [[1,1,0]]
                possible_ranges = [v for (k,v) in evt_e if k == original_range[1]][1] # all the possible ranges after the extinction from evt_e dict
                lost_range = [k for (k,v) in ranges_inv if v==next_event[1][2]][1] # extinction occurs in this region next_event[1][2], causing the original range to lose this region
                new_range = [v for v in possible_ranges if v != lost_range] # the new range for this branch after the event 
                new_label = [v for (k,v) in ranges_inv if k==new_range[1]][1] # the new label e.g. 1 for range [1,0,0] or 7 for [1,1,1]
                #
                # println(history_df)
                # assign the new range label to that branch idx (the branch is still in pool for either speciation or extinction)
                history_df[[row_chosen],:final_state]          .= new_label 
                # record the anagenetic change on that branch
                push!(anagen_df,(history_df[[row_chosen],:idx][1],curr_time,next_state[1],new_label))
                #
                # println("next event is extinction")
                # println(history_df)
                # println(anagen_df)
                # println(next_event[1][2])
                # xxx
            elseif next_event[1][1] == "e" && next_state[1] <= num_regions # if the next event is extinction on a branch with endemic range (Sep 5: already checked)
                # update the chosen branch's status 
                history_df[[row_chosen],:time_end]       .= curr_time
                history_df[[row_chosen],:extinct_status] .= 1
                history_df[[row_chosen],:split_status]   .= 0
                history_df[[row_chosen],:branch_length]  .= curr_time - deepcopy(history_df[[row_chosen],:time_start][1])
                history_df[[row_chosen],:event_type]     .= "e"
                # println(history_df)
                # xxx
            elseif next_event[1][1] == "d" # if the next event is dispersal (Sep 6: already checked)
                # update the chosen branch's status 
                original_range  = [k for (k,v) in ranges_inv if v==next_state[1]] # original range before dispersal event 
                possible_ranges =  [v for (k,v) in evt_d if k==original_range[1]][1] # all possible ranges following a dispersal event in that original_range
                new_range = [] # initialize empty range
                for i in 1:length(possible_ranges)
                    if possible_ranges[i][next_event[1][3]] == 1 # find the new range following a dispersal to a new region next_event[1][3] 
                        push!(new_range,possible_ranges[i])
                    end
                end
                new_label = [v for (k,v) in ranges_inv if k==new_range[1]][1] # label for the new range e.g. 5 for [1,0,1]
                # assign the new range label to that branch idx (the branch is still in pool for either speciation or extinction)
                history_df[[row_chosen],:final_state]          .= new_label 
                # record the anagenetic change on that branch
                push!(anagen_df,(history_df[[row_chosen],:idx][1],curr_time,next_state[1],new_label))
                #
                #  println("this is a dispersal event")
                #  println("current state = ", next_state[1])
                #  println("next event = ", next_event[1])
                #  println("original range = ", original_range)
                #  println("possible range = ", possible_ranges)
                #  println("new range = ",new_range)
                #  println("new label = ", new_label)
                #  println("history = ",history_df)
                #  println("anagenetic = ",anagen_df)
                #  xxx
            elseif next_event[1][1] == "b" # if the next event is between-region speciation (Sep 6: already checked)
                # println("history before = ",history_df)
                max_idx = maximum(history_df.idx) # find the maximum idx from current recoded branch labels in history_df
                push!(history_df,(max_idx+1,curr_time,missing,next_event[1][2],missing,missing,missing,missing,missing,missing)) # left range for child 1
                push!(history_df,(max_idx+2,curr_time,missing,next_event[1][3],missing,missing,missing,missing,missing,missing)) # right range for child 2
                # update the parent branch's status
                # P.S. row_chosen must be transformed into a vector instead of scalar (by doing [row_chosen])
                history_df[[row_chosen],:time_end]       .= curr_time
                history_df[[row_chosen],:extinct_status] .= 0
                history_df[[row_chosen],:split_status]   .= 1
                history_df[[row_chosen],:branch_length]  .= curr_time - deepcopy(history_df[[row_chosen],:time_start][1])
                history_df[[row_chosen],:child_1_idx]    .= deepcopy(max_idx+1)
                history_df[[row_chosen],:child_2_idx]    .= deepcopy(max_idx+2)
                history_df[[row_chosen],:event_type]     .= "b"
                #
                # println("this is between-region speciation")
                # println("event type = ",next_event[1] )
                # println("parent state = ",next_state[1])
                # println("history after = ",history_df)
                # xxx
            end
            
            ############################
            # update species count across regions in continuous time 
            ############################
            chosen_time = findall(x -> x <= curr_time, [v[1] for v in species_cts_time_region])
                for i in 1:num_regions
                    chosen_time = intersect(chosen_time,findall(x -> ismissing(x), [v[i+1] for v in species_cts_time_region]))
                end 
                 
                for j in 1:length(chosen_time)
                    global species_cts_time_region[chosen_time[j]][2:(num_regions+1)] = deepcopy(num_species_per_region)
                end
            
            ############################
            # update species count across ranges in continuous time 
            ############################
            chosen_time = findall(x -> x <= curr_time, [v[1] for v in species_cts_time_range])
                for i in 1:num_ranges
                    chosen_time = intersect(chosen_time,findall(x -> ismissing(x), [v[i+1] for v in species_cts_time_range]))
                end 
                                  
                for j in 1:length(chosen_time)
                    global species_cts_time_range[chosen_time[j]][2:(num_ranges+1)] = deepcopy(num_species_per_range)
                end

            ############################
            # # update num_species_per_region and num_species_per_range following an event (Sep 6: already checked)
            ############################
            # println("num_species before = ",num_species_per_range)
            # println("num_species region before = ",num_species_per_region)
            # update num_species_per_range following the previous event 
            for i in 1:num_ranges
                count_row = count(row -> row.final_state .==i && ismissing(row.extinct_status) && ismissing(row.split_status), eachrow(history_df))
                global num_species_per_range[i] = copy(count_row)
            end
            # update num_species_per_region following the previous event 
            new_species_count = fill(0, num_regions)
            for i in 1:length(num_species_per_range) # loop over each range in current num_species_per_range vector
                if num_species_per_range[i] != 0 # if there is a species (or more) with that range
                    range_key = [k for (k,v) in ranges_inv if v==i][1] # find the range key correspon to that range label
                    for j in 1:length(range_key) # loop over each region for that given range 
                        if range_key[j] !=0 # if there is a species (or more) in that region
                            global new_species_count[j] += num_species_per_range[i]*1 # record those species in that region
                        end
                    end
                end
            end

            global num_species_per_region = copy(new_species_count) # update species count in each region for the next event iteration
            # record num_species_per_region and num_species_per_range for this time step
            push!(range_count_history,deepcopy(num_species_per_range))
            range_count_history[end] = insert!(Float64[x for x in range_count_history[end]] ,1,curr_time)
            
            #
            push!(region_count_history,deepcopy(num_species_per_region))
            region_count_history[end] = insert!(Float64[x for x in region_count_history[end]] ,1,curr_time)
            #
            # update the number of currently alive lineages 
            global num_alive_lineage    = sum(num_species_per_range)
            # println("num alive species ",num_alive_lineage)
            # println("current time ",curr_time )
            # println("num species per range ",num_species_per_range)
            #
            # println("current history = ",history_df)
            # println("anagen history = ",anagen_df)
            # println("num_species after = ",num_species_per_range)
            # println("num_species region after = ",num_species_per_region)
            # println("history region count = ", region_count_history)
            # println("history range count = ", range_count_history)
            # println("num alive lineage = ",num_alive_lineage)
            # xxx 

        end
    #################################
    # # check if we stop simulation
    #################################
    
        ### Non tree extinction case
        if stop_condition == "time" && curr_time >= max_time && any(x -> x != 0, num_species_per_range)  # if condition on max_time and we don't have tree extinction
            global do_simulation = false
        elseif stop_condition == "numtaxa" && num_alive_lineage >= max_taxa && any(x -> x != 0, num_species_per_range)  # if condition on max_taxa and we don't have tree extinction
            global do_simulation = false
        ### Tree extinction case
        elseif all(x -> x == 0, num_species_per_range) && resimulate == false # we have tree extinction and resimulate == false
            global do_simulation = false 
        elseif all(x -> x == 0, num_species_per_range) && resimulate == true # we have tree extinction and  resimulate == true 
            global do_simulation = true
        end

    #################################
    # # post hoc tree treatment 
    #################################
    
        # stop simulation if condition met 
        if do_simulation == false && any(x -> x != 0, num_species_per_range) && stop_condition == "time" # if we don't have tree extinction and condition = "time"
            # # extend the branch length of all currently living lineages to the max_time 
            which_lineages = filter(row -> ismissing(row.extinct_status) && ismissing(row.split_status),history_df).idx
            # println(history_df)
            # println(which_lineages)
            for i in 1:length(which_lineages)
                history_df.time_end[which_lineages[i]]      = copy(max_time) # extend their external branches to max time
                # println(deepcopy(history_df.time_end[which_lineages[i]][1]))
                # println(deepcopy(history_df.time_start[which_lineages[i]][1]))
                # xxx
                history_df.branch_length[which_lineages[i]] = max_time - deepcopy(history_df.time_start[which_lineages[i]][1]) # compute the length of external branch
            end
            # println(history_df)
            # xxx
            # 
            # Fill the remaining rows with missing value of species_cts_time_region with num_spec_per_region from previous entry (since we just extend branch lengths here)
            missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_region])
            for i in 1:length(missing_rows)
                global species_cts_time_region[missing_rows[i]][2:(num_regions+1)] = deepcopy(num_species_per_region)
            end
            # Fill the remaining rows with missing value of  species_cts_time_range with num_spec_per_range from previous entry (since we just extend branch lengths here)
            missing_rows = findall(x -> ismissing(x), [v[2] for v in species_cts_time_range])
            for i in 1:length(missing_rows)
                global species_cts_time_range[missing_rows[i]][2:(num_ranges+1)] = deepcopy(num_species_per_range)
            end
        elseif do_simulation == false && any(x -> x != 0, num_species_per_range) && stop_condition == "numtaxa" # if we don't have tree extinction and condition = "numtaxa"
            break
        elseif do_simulation == false && all(x -> x == 0, num_species_per_range) # if we have tree extinction and resimulate == false
            println("tree extinction")
        elseif do_simulation == true && all(x -> x == 0, num_species_per_range)  # if we have tree extinction and resimulate == true
            println("tree extinction, " * "resimulate tree " * sim_id)
            # restart everything
            global num_species_per_range  = fill(0, num_ranges); # Empty vector to store number of species in each range
            global num_species_per_region = fill(0, num_regions); # Empty vector to store number of species in each region
            global start_time                         = 0.0; # Start time of the simulation (0 = past time)
            # global start_range                        = floor(Int64,rand(Uniform(1,num_ranges))); # Starting range state at root node e.g. 1,2,...,7 for 3 region
            global start_range                        = 1;
            global num_alive_lineage                  = 1; # Start with one species (stem branch)
            num_species_per_range[start_range] = 1;
            #
            global range_key_root                          = [k for (k,v) in ranges_inv if v==start_range]; #range key that correspond to range label at the root e.g. (1,1,0) for label 4
            global reg_idx_root                            = findall(range_key_root[1] .==1 ); # find the region indexes that correspond to that range at the root e.g. region_ind = 1,2 for range (1,1,0)
            num_species_per_region[reg_idx_root]    .= 1; # update the number of species currently alive in each region e.g. 1 species with range {A,B} means 1 species in A and 1 species in B
            #
            global range_count_history  = [] # record num_species_per_range across time  
            global region_count_history = [] # record num_species_per_region across time 
            #
            push!(range_count_history,deepcopy(num_species_per_range))
            range_count_history[1] = insert!(Float64[x for x in range_count_history[1]] ,1,start_time)
            #
            push!(region_count_history,deepcopy(num_species_per_region))
            region_count_history[1] = insert!(Float64[x for x in region_count_history[1]] ,1,start_time)
            #
            global species_cts_time_region = [] # initialize array of vectors containing information about num_species across regions in cts time
            for i in 1:timebin
                global added_vec     = fill(missing,num_regions)
                global added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
                insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
                #
                push!(species_cts_time_region,added_vec_mod)
            end 
            species_cts_time_region[1][2:(num_regions+1)]= deepcopy(num_species_per_region)
            #
            global species_cts_time_range = [] # initialize array of vectors containing information about num_species across ranges in cts time
            for i in 1:timebin
                global added_vec     = fill(missing,num_ranges)
                global added_vec_mod = Union{Any,Missing}[added_vec...] # allow another value to be added to a vector that contains missing type
                insert!(added_vec_mod,1,cts_time_pts[i]) # insert each continuous time point to each vector
                #
                push!(species_cts_time_range,added_vec_mod)
            end 
            species_cts_time_range[1][2:(num_ranges+1)]=deepcopy(num_species_per_range)        
            #
            curr_time = start_time; # Initialize current time of simulation p.s. time starts from 0 (past)
            global history_df = DataFrame(idx = [],time_start = [],time_end = [],final_state = [],extinct_status = [],split_status = [],branch_length = [],child_1_idx = [],child_2_idx = [],event_type = []) 
            push!(history_df,(1,start_time,missing,start_range,missing,missing,missing,missing,missing,missing)) # add first event at the root node (root node)
            global anagen_df = DataFrame(idx = [],time = [], from_state = [], to_state = [])
        end
        
    end

    ## Enable this if want to validate the simulation
    # println("final history = ",history_df)
    # println("final anagenetic history = ", anagen_df)
    # println("region count over time = ",region_count_history)
    # println("range count over time",range_count_history)
    # println("final num alive lineage = ",num_alive_lineage)
    # # xxx
    # println("length history_df = ",nrow(history_df))

    #########################
    # Character data at tips 
    #########################
    tip_range_df      = filter([:extinct_status, :split_status] => (l, v) -> ismissing(l) && ismissing(v), history_df) # get tip state information
    tip_character_dat = Matrix{Any}(undef, nrow(tip_range_df), num_regions+1) # empty matrix to store tip state and tip label 
    for i in 1:nrow(tip_range_df)
        region_arr = deepcopy([k for (k,v) in ranges_inv if v==tip_range_df.final_state[i]][1]) # region presence/absence that correspond to range in range_inv
        tip_character_dat[i,1] = tip_range_df.idx[i] # assign tip index 
        for j in 1:num_regions
            tip_character_dat[i,(j+1)] = region_arr[j] # assign region presence/absence
        end
    end
    tip_character_df = DataFrame(tip_character_dat,:auto) #convert to dataframe

    #########################
    # Rate labels and other summary statistics
    #########################
    param_label_df =  DataFrame(base_w = rho_w,base_e = rho_e,base_d = rho_d,base_b = rho_b,scalar_w = w_d,
                                scalar_e = e_d, scalar_d_i = d_d_home, scalar_d_j = d_d_away, scalar_b = b_d, num_alive_spec = num_alive_lineage);
    
    #########################
    # Write Newick string 
    #########################
    newick_str   = "" # initialize an empty string
    #
    # starting from the root branch 
    if ismissing(history_df.split_status[1]) # if the root branch survive till present without any speciation
        global newick_str = "(sp." * string(history_df.idx[1]) * ":" * string(history_df.branch_length[1])*")"
    else
        if history_df.split_status[1] == 1 # if the root branch ends with a speciation event, do (branch2., branch3.):length_1
            global newick_str = "(branch2.,branch3.):" * string(history_df.branch_length[1])
        elseif history_df.split_status[1] == 0 # if the root branch ends with a extinction event
            global newick_str = "(sp." * string(history_df.idx[1]) * ":" * string(history_df.branch_length[1])*")"
        end
    end

    # loop over branches descending from the root branch
    for i in 2:nrow(history_df)
        if ismissing(history_df.split_status[i]) # if branch i is a tip branch
            prev_label  = "branch"*string(i)*"."
            next_label  = "sp." * string(i) * ":" * string(history_df.branch_length[i])
            global newick_str = replace(newick_str, prev_label => next_label)
        else # if branch i is not a tip branch
            if history_df.split_status[i] == 1 # if branch i ends with a speciation event, replace string i with its child indexes j,k e.g. "branchi." => "(branchj.,branchk.):length_i"
                prev_label = "branch"*string(i)*"."
                next_label = "(" * "branch"*string(history_df.child_1_idx[i]) * "." * "," * "branch"*string(history_df.child_2_idx[i]) * "." * "):" * string(history_df.branch_length[i])
                global newick_str = replace(newick_str,prev_label => next_label)
            elseif history_df.split_status[i] == 0 # if branch i ends with an extinction event
                prev_label = "branch"*string(i)*"."
                next_label = "sp." * string(i) * ":" * string(history_df.branch_length[i])
                global newick_str = replace(newick_str,prev_label => next_label)
            end
        end 
    end
    global newick_str = newick_str * ";"
    # println("newick = ",newick_str)
    # xxx
    # end

    ################## 
    # # Write outputs
    #################
    
    # # write tree as a newick string
    # println("newick = ",newick_str)

    tree_file = "complete_tree_" * sim_id * ".tre"
    println("done simulating tree "*  sim_id)
    touch(tree_file) 
    open(tree_file,"w") do file
        write(file,newick_str)
    end
    
    # # write tip data to file (as csv)
    CSV.write("tip_data_tree_"*sim_id*".dat.csv",tip_character_df;header=false)
    # println(history_df)
    # println(tip_range_df)
    # println(tip_character_dat)
    # # write parameters to file

    # # convert each summary statistics to data frame 
    dat_reg_cts   = DataFrame(time = [x[1] for x in species_cts_time_region],region_1 = [x[2] for x in species_cts_time_region],
                            region_2 = [x[3] for x in species_cts_time_region],region_3 = [x[4] for x in species_cts_time_region])

    dat_range_cts = DataFrame(time = [x[1] for x in species_cts_time_range],range_1 = [x[2] for x in species_cts_time_range],
                            range_2 = [x[3] for x in species_cts_time_range],range_3 = [x[4] for x in species_cts_time_range],
                            range_4 = [x[5] for x in species_cts_time_range],range_5 = [x[6] for x in species_cts_time_range],
                            range_6 = [x[7] for x in species_cts_time_range],range_7 = [x[8] for x in species_cts_time_range])

    dat_reg   = DataFrame(time = [x[1] for x in region_count_history],region_1 = [x[2] for x in region_count_history],
                            region_2 = [x[3] for x in region_count_history],region_3 = [x[4] for x in region_count_history])

    dat_range = DataFrame(time = [x[1] for x in range_count_history],range_1 = [x[2] for x in range_count_history],
                            range_2 = [x[3] for x in range_count_history],range_3 = [x[4] for x in range_count_history],
                            range_4 = [x[5] for x in range_count_history],range_5 = [x[6] for x in range_count_history],
                            range_6 = [x[7] for x in range_count_history],range_7 = [x[8] for x in range_count_history])   

    # # write summary statistics and rates to file
    CSV.write("param_label_tree_"*sim_id*".labels.csv",param_label_df) # write rate parameters and number of alive tips
    CSV.write("anagenetic_changes_tree_"*sim_id*".csv",anagen_df) # write anagenetic history along branches over event times 
    CSV.write("species_count_region_tree_"*sim_id*".csv",dat_reg) # write species count across regions over event times
    CSV.write("species_count_range_tree_"*sim_id*".csv",dat_range) # write species count across ranges over event times
    CSV.write("cts_species_count_region_tree_"*sim_id*".csv",dat_reg_cts) # write species count across regions over cts times
    CSV.write("cts_species_count_range_tree_"*sim_id*".csv",dat_range_cts) # write species count across ranges over cts times
    CSV.write("cladogenetic event history_"*sim_id*".csv",history_df) # write cladogenetic event history over event times 
# end