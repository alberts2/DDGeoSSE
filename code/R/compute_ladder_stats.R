# this script is similar to the original plot_ddgeosse_treestats.R 
# but we are computing corrected Colless, rather than beta.

# 8 June
# use the y-axis limit for the corrected colless plot 
library(ape)
library(treestats)
library(geiger)
library(gridExtra)
library(grid)
library(ggplot2)

which_dat = "small" # small or larger dd interval
fp      = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
tree_fp = paste0(fp,"data/tree_stats/smaller_interval/data_treestats_small/")
# tree_fp = paste0(fp,"data/tree_stats/larger_interval/data_treestats_big/")
dat_fp  = paste0(fp,"data/tree_stats/smaller_interval/")
# dat_fp  = paste0(fp,"data/tree_stats/larger_interval/")
out_fp    = paste0(fp,"plot/sim_treestats/")

num_sim = 50 
mesh = read.table(paste0(dat_fp,'mesh_tree_stats_50_small.csv'),sep = ",",header = T)
# mesh = read.table(paste0(dat_fp,'mesh_tree_stats_50_big.csv'),sep = ",",header = T)

mesh$idx = seq(1,nrow(mesh))

# wd
mesh_wd              = mesh[mesh$w_d != 0 & mesh$e_d == 0 & mesh$d_d_away == 0,]
mesh_wd$avg_ladder   = rep(NA,nrow(mesh_wd))

for (i in 1:nrow(mesh_wd)){
  ladder_vec = c() 
  for (j in (((mesh_wd$idx[i]-1)*num_sim)+1):(mesh_wd$idx[i]*num_sim)){
    tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
    tre_extant = drop.extinct(tre)
    ladder_vec = append(ladder_vec,avg_ladder(tre_extant))
    print(paste0("finish tree ", j))
  }
  mesh_wd$avg_ladder[i] = mean(ladder_vec,na.rm = TRUE)
  print(paste0("finish i = ",i))
}

if (which_dat == "small"){
  mesh_wd$colless_corr   = rep(NA,nrow(mesh_wd))
  
  for (i in 1:nrow(mesh_wd)){
    colless_vec = c() 
    for (j in (((mesh_wd$idx[i]-1)*num_sim)+1):(mesh_wd$idx[i]*num_sim)){
      tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
      tre_extant = drop.extinct(tre)
      colless_vec = append(colless_vec,colless_corr(tre_extant))
      print(paste0("finish tree ", j))
    }
    mesh_wd$colless_corr[i] = mean(colless_vec,na.rm = TRUE)
    print(paste0("finish i = ",i))
  }
}

# ed
mesh_ed              = mesh[mesh$w_d == 0 & mesh$e_d != 0 & mesh$d_d_away == 0,]
mesh_ed$avg_ladder   = rep(NA,nrow(mesh_ed))

for (i in 1:nrow(mesh_ed)){
  ladder_vec = c() 
  for (j in (((mesh_ed$idx[i]-1)*num_sim)+1):(mesh_ed$idx[i]*num_sim)){
    tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
    char_dat =  read.table(paste0(tree_fp,'sim.',j,".extant.dat.csv"),sep = ",",header = T)
    # ignore tree with only 1 extant taxa 
    if (nrow(char_dat) != 1){
      tre_extant = drop.extinct(tre)
      # Normalize colless does not work for tree with only 2 taxa, so ignore that tree too
      ladder_vec = append(ladder_vec,avg_ladder(tre_extant))
      print(paste0("finish tree ", j))
    }
  }
  mesh_ed$avg_ladder[i] = mean(ladder_vec,na.rm = TRUE)
  print(paste0("finish i = ",i))
}

if (which_dat == "small"){
  mesh_ed$colless_corr   = rep(NA,nrow(mesh_ed))
  
  for (i in 1:nrow(mesh_ed)){
    colless_vec = c() 
    for (j in (((mesh_ed$idx[i]-1)*num_sim)+1):(mesh_ed$idx[i]*num_sim)){
      tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
      tre_extant = drop.extinct(tre)
      colless_vec = append(colless_vec,colless_corr(tre_extant))
      print(paste0("finish tree ", j))
    }
    mesh_ed$colless_corr[i] = mean(colless_vec,na.rm = TRUE)
    print(paste0("finish i = ",i))
  }
}

# dd
mesh_dd              = mesh[mesh$w_d == 0 & mesh$e_d == 0 & mesh$d_d_away != 0,]
mesh_dd$avg_ladder   = rep(NA,nrow(mesh_dd))

for (i in 1:nrow(mesh_dd)){
  ladder_vec = c() 
  for (j in (((mesh_dd$idx[i]-1)*num_sim)+1):(mesh_dd$idx[i]*num_sim)){
    tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
    char_dat =  read.table(paste0(tree_fp,'sim.',j,".extant.dat.csv"),sep = ",",header = T)
    # ignore tree with only 1 extant taxa 
    if (nrow(char_dat) != 1){
      tre_extant = drop.extinct(tre)
      ladder_vec = append(ladder_vec,avg_ladder(tre_extant))
    }
  }
  mesh_dd$avg_ladder[i] =mean(ladder_vec,na.rm = TRUE)
}

if (which_dat == "small"){
  mesh_dd$colless_corr   = rep(NA,nrow(mesh_dd))
  
  for (i in 1:nrow(mesh_dd)){
    colless_vec = c() 
    for (j in (((mesh_dd$idx[i]-1)*num_sim)+1):(mesh_dd$idx[i]*num_sim)){
      tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
      tre_extant = drop.extinct(tre)
      colless_vec = append(colless_vec,colless_corr(tre_extant))
      print(paste0("finish tree ", j))
    }
    mesh_dd$colless_corr[i] = mean(colless_vec,na.rm = TRUE)
    print(paste0("finish i = ",i))
  }
}

# save file
write.csv(mesh_wd,file = paste0(out_fp, "ladder_wd_small.csv"),row.names = FALSE)

write.csv(mesh_ed,file = paste0(out_fp, "ladder_ed_small.csv"),row.names = FALSE)

write.csv(mesh_dd,file = paste0(out_fp, "ladder_dd_small.csv"),row.names = FALSE)

# write.csv(mesh_wd,file = paste0(out_fp, "ladder_wd_big.csv"),row.names = FALSE)
# 
# write.csv(mesh_ed,file = paste0(out_fp, "ladder_ed_big.csv"),row.names = FALSE)
# 
# write.csv(mesh_dd,file = paste0(out_fp, "ladder_dd_big.csv"),row.names = FALSE)

