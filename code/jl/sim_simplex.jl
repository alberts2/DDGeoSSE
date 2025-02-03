#####################
# Set file system
#####################
fp       = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
input_fp = fp*"code/jl/"
out_fp   = fp*"data/tree_stats/"
#####################
# Load julia packages
#####################
using RCall # Package to call R function
using DataFrames # Dataframes package
using CSV # Package to write to csv

#####################
# Load R packages
#####################
R"
options(rlibrary_normalize=FALSE)
library(treestats) # Package to compute treestats
library(ape) # Package to read tree file
library(geiger) # Package to produce extant-only tree
"

#####################
# Create simplex grid
#####################
R"
num_steps = 10 # number of each individual point on mesh grid
# num_steps = 5 # number of each individual point on mesh grid
w_d = d_d_away <- round(seq(-2,0,length.out = num_steps),2)
e_d            <- round(seq(0,2,length.out = num_steps),2)
# w_d = d_d_away <- round(seq(-0.1,0,length.out = num_steps),2)
# e_d            <- round(seq(0,0.1,length.out = num_steps),2)
mesh           <- expand.grid(w_d,d_d_away,e_d)
mesh           <- cbind(mesh,rep(1,nrow(mesh)),rep(1,nrow(mesh)),rep(0,nrow(mesh)),rep(0,nrow(mesh)),rep(0,nrow(mesh)),rep(0,nrow(mesh))) #placeholder for beta, gamma, mean range size, var range size,#discarded trees columns, mean extant taxa
colnames(mesh) <- c('w_d','d_d_away','e_d','beta','gamma','mean_range_size','var_range_size','prop_invalid_trees','mean_extant_taxa') 
mesh           <- unique(mesh[,1:9]) # remove duplicated rows from rounding
mesh <- mesh[!(mesh$w_d == 0 & mesh$d_d_away == 0 & mesh$e_d == 0), ] # remove entry that correspond to standard geosse
"

#####################
# Simulate batch of trees for each point vector in mesh grid
#####################
mesh = @rget mesh # pass mesh grid created in r to julia 

# Define the arguments
num_sim     = 50 # number of trees per combination
num_regions = 3 # num of regions, matched with setting in sim_one_variable.jl
num_row     = nrow(mesh) # number of parameter combinations 
arg1        = "simulate"
# arg1        = "resimulate"
arg2        = "sim"


for i in 1:num_row #loop over each parameter combination on mesh grid 
    arg3 = (i-1)*num_sim + 1 #starting tree index
    arg4 = num_sim # number of trees to simulate
    arg5 = Float64(mesh[i,1]) # w_d
    arg6 = Float64(mesh[i,2]) # d_d_away
    arg7 = Float64(mesh[i,3]) # e_d
    # Construct the command
    script_path = "sim_one_variable.jl"
    cmd = `julia $script_path $arg1 $arg2 $arg3 $arg4 $arg5 $arg6 $arg7`
    # Run sim_one_func.jl
    run(cmd)
    # Get beta stat for each tree in the tree set for that parameter combination
    @rput input_fp # pass input_fp from julia to R
    @rput num_sim  # pass num_sim from julia to R
    @rput i        # pass i from julia to R
    @rput num_regions # pass num_regions from julia to R
    R"
    beta_stat           = c()
    gamma_stat          = c()
    num_extant_taxa     = c()
    mean_range_size_vec = c()
    var_range_size_vec  = c()
    count               = 0
    for (j in ((i-1)*num_sim + 1):(i*num_sim)){
            tree_dat           = read.tree(paste0(input_fp,'simulate/sim.',j,'.tre'))
            tree_dat$root.edge = 0 # set root edge length equal 0
            extant_dat         = try(drop.extinct(tree_dat), silent = TRUE)
            #
            range_dat = read.table(paste0(input_fp,'simulate/sim.',j,'.extant.dat.csv'),sep = ',',header = T)
            if(inherits(extant_dat, 'try-error')){ # if complete tree only has 1 tip (either alive or ext), then dont update extant tree and dont compute tree stats
                extant_dat = NA
            }
            # Compute tree stats
            if (all(!is.na(extant_dat)) == TRUE && (extant_dat$Nnode > 2)){
                beta_val        = beta_statistic(extant_dat)
                gamma_val       = gamma_statistic(extant_dat)
                mean_range_size = mean(rowSums(range_dat[,2:(num_regions+1)] == 1)) # compute average range size across species for tree j
                var_range_size  = var(rowSums(range_dat[,2:(num_regions+1)] == 1))  # compute variance range size across species for tree j
                #
                beta_stat           = append(beta_stat,beta_val)
                gamma_stat          = append(gamma_stat,gamma_val)
                num_extant_taxa     = append(num_extant_taxa,length(extant_dat$tip.label))
                mean_range_size_vec = append(mean_range_size_vec,mean_range_size)
                var_range_size_vec  = append(var_range_size_vec,var_range_size)

            } else { # if we have invalid tree for computing beta
                beta_stat  = append(beta_stat,NA)
                gamma_stat = append(gamma_stat,NA)
                count      = count + 1 
            }
    }
    aver_beta        = mean(beta_stat,na.rm = T)
    aver_gamma       = mean(gamma_stat,na.rm = T)
    mean_extant_taxa = mean(num_extant_taxa,na.rm = T)
    prop_invalid     = count / num_sim
    mean_range_size  = mean(mean_range_size_vec,na.rm = T) # compute average range size across trees from combination i 
    var_range_size   = mean(var_range_size_vec,na.rm = T)  # compute variance range size across trees from combination i 
    "
    aver_beta        = @rget aver_beta
    aver_gamma       = @rget aver_gamma
    prop_invalid     = @rget prop_invalid
    mean_range_size  = @rget mean_range_size
    var_range_size   = @rget var_range_size
    mean_extant_taxa = @rget mean_extant_taxa
    #
    mesh[i,4] = aver_beta
    mesh[i,5] = aver_gamma
    mesh[i,6] = mean_range_size
    mesh[i,7] = var_range_size
    mesh[i,8] = prop_invalid
    mesh[i,9] = mean_extant_taxa
end

#####################
# Save output
#####################
CSV.write(joinpath(out_fp,"mesh_tree_stats_"*string(num_sim)*".csv"),mesh)

#####################
# Load output
#####################



