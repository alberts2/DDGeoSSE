# This script simulates standard geosse trees and computes some tree statistics from them

library(treestats)

##### FILE SYSTEM #####
#######################
fp = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
out_fp = paste0(fp,"data/sim_geosse/diversitree/") #output 

##### MODEL PARAMETERS #####
#######################
# w_A, w_B, b_A_B, e_A, e_B, d_AB, d_BA
pars <- c(1.0, 0.5, 0.5, 0.2, 0.4, 2.5, 0.5)

num_sim           = 1000
max_time          = 7
root_state        = 0

balance_diversitree_df = data.frame(colless = NA, colless_corr = NA, colless_quad = NA,
                                    sackin = NA, beta = NA, stairs = NA, stairs2 = NA,
                                    blum_norm = NA, rogers_norm = NA, pitchforks_norm = NA
                                    )

# Note: crown age is the starting time below the first split (i.e. assumes that root branch length = 0)
length_diversitree_df = data.frame(crown_age = NA, gamma = NA, pigot_rho = NA)


for (i in 1:num_sim){
  # simulate a single phylogeny until max_time unit time and root state x0 (region A)
  phy <- tree.geosse(pars, max.t=max_time, x0=root_state)
  while (length(phy$tip.label) == 0 || length(phy$tip.label) == 2){ # if either we have tree ext or tree with only 2 extant taxa (cant compute beta)
    phy <- tree.geosse(pars, max.t=max_time, x0=root_state) # then re-simulate tree until valid
  }
    # save the phylogeny and information
    save(phy, file = paste0(out_fp,"diversitree_phylo_dat_",i,".RData"))
    write.tree(phy, file = paste0(out_fp,"diversitree_phylo_",i,".newick"))
    # compute the balance statistics and store them 
    balance_diversitree_df[i,] = c(colless(phy),colless_corr(phy),colless_quad(phy),
                                   sackin(phy),beta_statistic(phy),stairs(phy),stairs2(phy),
                                   blum(phy,normalization = T),rogers(phy,normalization = "tips"),
                                   pitchforks(phy,normalization = "tips"))
    # compute the branch length statistics and store them
    length_diversitree_df[i,] = c(crown_age(phy),treestats::gamma_statistic(phy),pigot_rho(phy))
}

write.table(balance_diversitree_df,file = paste0(out_fp,"diversitree_dat",".csv"), sep = ",", col.names = T,row.names = F)
write.table( length_diversitree_df,file = paste0(out_fp,"diversitree_length_dat",".csv"), sep = ",", col.names = T,row.names = F)

