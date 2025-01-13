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
w_d = d_d_away <- round(seq(-2,0,length.out = num_steps),2)
e_d            <- round(seq(0,2,length.out = num_steps),2)
mesh           <- expand.grid(w_d,d_d_away,e_d)
mesh           <- cbind(mesh,rep(1,nrow(mesh)),rep(1,nrow(mesh)))
colnames(mesh) <- c('w_d','d_d_away','e_d','beta','gamma') 
mesh           <- unique(mesh[,1:5]) # remove duplicated rows from rounding
mesh <- mesh[!(mesh$w_d == 0 & mesh$d_d_away == 0 & mesh$e_d == 0), ] # remove entry that correspond to standard geosse
"

#####################
# Simulate batch of trees for each point vector in mesh grid
#####################
mesh = @rget mesh # pass mesh grid created in r to julia 

# Define the arguments
num_sim = 10 # number of trees per combination
num_row = nrow(mesh) # number of parameter combinations 
arg1 = "simulate"
arg2 = "sim"


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
    R"
    beta_stat  = c()
    gamma_stat = c()
    for (j in ((i-1)*num_sim + 1):(i*num_sim)){
            print(j)
            tree_dat = read.tree(paste0(input_fp,'simulate/sim.',j,'.tre'))
            tree_dat$root.edge = 0 # set root edge length equal 0
            extant_dat = try(drop.extinct(tree_dat), silent = TRUE)
            if(inherits(extant_dat, 'try-error')){ # if complete tree only has 1 tip (either alive or ext), then dont update extant tree and dont compute tree stats
                extant_dat = NA
            }
            # Compute tree stats
            if (all(!is.na(extant_dat)) == TRUE && (extant_dat$Nnode > 2)){
                beta_val = beta_statistic(extant_dat)
                gamma_val = gamma_statistic(extant_dat)
                beta_stat = append(beta_stat,beta_val)
                gamma_stat = append(gamma_stat,gamma_val)
            } else { # if we have invalid tree for computing beta
                beta_stat = append(beta_stat,NA)
                gamma_stat = append(gamma_stat,NA)
            }
    }
    aver_beta = mean(beta_stat,na.rm = T)
    aver_gamma = mean(gamma_stat,na.rm = T)
    "
    aver_beta = @rget aver_beta
    aver_gamma = @rget aver_gamma
    #
    mesh[i,4] = aver_beta
    mesh[i,5] = aver_gamma
end

#####################
# Save output
#####################
CSV.write(joinpath(out_fp,"mesh_tree_stats_"*string(num_sim)*".csv"),mesh)



