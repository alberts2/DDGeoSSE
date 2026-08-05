# this script is similar to the original plot_ddgeosse_treestats.R 
# but we are computing corrected Colless, rather than beta.

# 8 June
# use the y-axis limit for the corrected colless plot 

# NOTE: tree with <= 3 tips does not apply for I statistic
# tree with <= 2 tips does not apply for beta statistic
library(ape)
library(treestats)
library(geiger)
library(gridExtra)
library(grid)
library(ggplot2)

fp      = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
tree_fp = paste0(fp,"data/tree_stats/smaller_interval/data_treestats_small/")
# tree_fp = paste0(fp,"data/tree_stats/larger_interval/data_treestats_big/")
dat_fp  = paste0(fp,"data/tree_stats/smaller_interval/")
# dat_fp  = paste0(fp,"data/tree_stats/larger_interval/")
out_fp          = paste0(fp,"plot/sim_treestats/")

num_sim = 50 
mesh = read.table(paste0(dat_fp,'mesh_tree_stats_50_small.csv'),sep = ",",header = T)
# mesh = read.table(paste0(dat_fp,'mesh_tree_stats_50_big.csv'),sep = ",",header = T)

mesh$idx = seq(1,nrow(mesh))

# wd
mesh_wd              = mesh[mesh$w_d != 0 & mesh$e_d == 0 & mesh$d_d_away == 0,]
# mesh_wd = mesh_wd[, c("idx", "w_d", "d_d_away", "e_d",
#                     "mean_range_size", "var_range_size",
#                     "prop_invalid_trees", "mean_extant_taxa",
#                     "log_mean_extant_taxa"
#                   )]
# for smaller case
mesh_wd = mesh_wd[, c("idx", "w_d", "d_d_away", "e_d",
                      "mean_range_size", "var_range_size",
                      "prop_invalid_trees", "mean_extant_taxa"
                  )]
dat_wd_stats = data.frame()

for (i in 1:nrow(mesh_wd)) {
  
  tmp = data.frame()
  
  for (j in (((mesh_wd$idx[i]-1)*num_sim + 1):(mesh_wd$idx[i]*num_sim))) {
    
    tre = read.tree(paste0(tree_fp, "sim.", j, ".tre"))
    tre_extant = drop.extinct(tre)
    # calculate treestats for that tre
    get_stats = calc_all_stats(tre_extant)
    # add to temporary data frame
    tmp = rbind(tmp, as.data.frame(as.list(get_stats)))
    print(paste0("finish tree ",j))
  }
  ## Mean over all simulations (j) for this i
  dat_wd_stats = rbind(dat_wd_stats,as.data.frame(as.list(colMeans(tmp))))
}
# combine the columns
mesh_wd = cbind(mesh_wd,dat_wd_stats)

# ed
mesh_ed              = mesh[mesh$w_d == 0 & mesh$e_d != 0 & mesh$d_d_away == 0,]
# mesh_ed = mesh_ed[, c("idx", "w_d", "d_d_away", "e_d",
#                       "mean_range_size", "var_range_size",
#                       "prop_invalid_trees", "mean_extant_taxa",
#                       "log_mean_extant_taxa"
#                   )]
# for smaller case 
mesh_ed = mesh_ed[, c("idx", "w_d", "d_d_away", "e_d",
                      "mean_range_size", "var_range_size",
                      "prop_invalid_trees", "mean_extant_taxa"
                  )]

dat_ed_stats <- data.frame()

for (i in seq_len(nrow(mesh_ed))) {
  tmp = data.frame()
  #
  j_start = (mesh_ed$idx[i] - 1) * num_sim + 1
  j_end   = mesh_ed$idx[i] * num_sim
  #
  for (j in j_start:j_end) {
    char_dat = read.table(paste0(tree_fp, "sim.", j, ".extant.dat.csv"),
                          sep = ",",header = TRUE
                          )
    # Ignore trees with only one extant taxon
    if (nrow(char_dat) > 1) {
      
      tre = read.tree(paste0(tree_fp, "sim.", j, ".tre"))
      # get extant tree
      tre_extant = drop.extinct(tre)
      # calculate tree stats for that tree
      get_stats = calc_all_stats(tre_extant)
      # add to dataframe
      tmp = rbind(tmp,as.data.frame(as.list(get_stats)))
    } else {
      print(paste0("Skipping tree ", j, ": only one extant taxon"))
    }
    print(paste0("Finished tree ", j))
  }
  # Calculate means only if at least one valid tree was retained
  if (nrow(tmp) > 0) {
    print("Computing means...")
    mean_stats = colMeans(tmp,na.rm = TRUE)
    #
    print("Appending row...")
    dat_ed_stats = rbind(dat_ed_stats,as.data.frame(as.list(mean_stats)))
  } else {
    warning(
      "No valid trees for mesh row ", i,
      " with idx = ", mesh_ed$idx[i]
    )
    # Add an NA row if no valid trees exist for this parameter combination
    if (ncol(dat_ed_stats) > 0) {
      dat_ed_stats[nrow(dat_ed_stats) + 1, ] <- NA_real_
    }
  }
  print(paste0("Finished i =", i))
}
# combine the columns
mesh_ed = cbind(mesh_ed,dat_ed_stats)

# dd
mesh_dd              = mesh[mesh$w_d == 0 & mesh$e_d == 0 & mesh$d_d_away != 0,]
# mesh_dd = mesh_dd[, c("idx", "w_d", "d_d_away", "e_d",
#                       "mean_range_size", "var_range_size",
#                       "prop_invalid_trees", "mean_extant_taxa",
#                       "log_mean_extant_taxa"
#                   )]
# for smaller case 
mesh_dd = mesh_dd[, c("idx", "w_d", "d_d_away", "e_d",
                      "mean_range_size", "var_range_size",
                      "prop_invalid_trees", "mean_extant_taxa"
                  )]

dat_dd_stats = data.frame()

for (i in seq_len(nrow(mesh_dd))) {
  tmp = data.frame()
  #
  j_start = (mesh_dd$idx[i] - 1) * num_sim + 1
  j_end   = mesh_dd$idx[i] * num_sim
  #
  for (j in j_start:j_end) {
    char_dat = read.table(paste0(tree_fp, "sim.", j, ".extant.dat.csv"),
                          sep = ",",header = TRUE
    )
    # Ignore trees with only one extant taxon
    if (nrow(char_dat) > 1) {
      
      tre = read.tree(paste0(tree_fp, "sim.", j, ".tre"))
      # get extant tree
      tre_extant = drop.extinct(tre)
      # calculate tree stats for that tree
      get_stats = calc_all_stats(tre_extant)
      # add to dataframe
      tmp = rbind(tmp,as.data.frame(as.list(get_stats)))
    } else {
      print(paste0("Skipping tree ", j, ": only one extant taxon"))
    }
    print(paste0("Finished tree ", j))
  }
  # Calculate means only if at least one valid tree was retained
  if (nrow(tmp) > 0) {
    print("Computing means...")
    mean_stats = colMeans(tmp,na.rm = TRUE)
    #
    print("Appending row...")
    dat_dd_stats = rbind(dat_dd_stats,as.data.frame(as.list(mean_stats)))
  } else {
    warning(
      "No valid trees for mesh row ", i,
      " with idx = ", mesh_dd$idx[i]
    )
    # Add an NA row if no valid trees exist for this parameter combination
    if (ncol(dat_dd_stats) > 0) {
      dat_dd_stats[nrow(dat_dd_stats) + 1, ] <- NA_real_
    }
  }
  print(paste0("Finished i =", i))
}
# combine the columns
mesh_dd = cbind(mesh_dd,dat_dd_stats)

# # save file
# write.csv(mesh_wd,file = paste0(out_fp, "mesh_wd_big.csv"),row.names = FALSE)
# 
# write.csv(mesh_ed,file = paste0(out_fp, "mesh_ed_big.csv"),row.names = FALSE)
# 
# write.csv(mesh_dd,file = paste0(out_fp, "mesh_dd_big.csv"),row.names = FALSE)

# save file
write.csv(mesh_wd,file = paste0(out_fp, "mesh_wd_small.csv"),row.names = FALSE)

write.csv(mesh_ed,file = paste0(out_fp, "mesh_ed_small.csv"),row.names = FALSE)

write.csv(mesh_dd,file = paste0(out_fp, "mesh_dd_small.csv"),row.names = FALSE)
