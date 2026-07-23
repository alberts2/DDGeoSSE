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

# append the corrected colless indices 
for (i in 1:nrow(mesh)){
  coll_vec = c() 
  for (j in (((mesh$idx[i]-1)*num_sim)+1):(mesh$idx[i]*num_sim)){
    tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
    char_dat =  read.table(paste0(tree_fp,'sim.',j,".extant.dat.csv"),sep = ",",header = T)
    # ignore tree with only 1 extant taxa 
    if (nrow(char_dat) != 1){
      tre_extant = drop.extinct(tre)
      # Normalize colless does not work for tree with only 2 taxa, so ignore that tree too
      coll_vec = append(coll_vec,colless_corr(tre_extant))
    }
  }
  mesh$colless_norm[i] = mean(coll_vec,na.rm = TRUE)
}

# wd
mesh_wd              = mesh[mesh$w_d != 0 & mesh$e_d == 0 & mesh$d_d_away == 0,]
mesh_wd$colless_norm = rep(NA,nrow(mesh_wd))

for (i in 1:nrow(mesh_wd)){
  coll_vec = c() 
  for (j in (((mesh_wd$idx[i]-1)*num_sim)+1):(mesh_wd$idx[i]*num_sim)){
    tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
    tre_extant = drop.extinct(tre)
    coll_vec = append(coll_vec,colless_corr(tre_extant))
  }
  mesh_wd$colless_norm[i] = mean(coll_vec,na.rm = TRUE)
}

# ed
mesh_ed              = mesh[mesh$w_d == 0 & mesh$e_d != 0 & mesh$d_d_away == 0,]
mesh_ed$colless_norm = rep(NA,nrow(mesh_ed))

for (i in 1:nrow(mesh_ed)){
  coll_vec = c() 
  for (j in (((mesh_ed$idx[i]-1)*num_sim)+1):(mesh_ed$idx[i]*num_sim)){
    tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
    char_dat =  read.table(paste0(tree_fp,'sim.',j,".extant.dat.csv"),sep = ",",header = T)
    # ignore tree with only 1 extant taxa 
    if (nrow(char_dat) != 1){
      tre_extant = drop.extinct(tre)
      # Normalize colless does not work for tree with only 2 taxa, so ignore that tree too
      coll_vec = append(coll_vec,colless_corr(tre_extant))
    }
  }
  mesh_ed$colless_norm[i] = mean(coll_vec,na.rm = TRUE)
}

# dd
mesh_dd              = mesh[mesh$w_d == 0 & mesh$e_d == 0 & mesh$d_d_away != 0,]
mesh_dd$colless_norm = rep(NA,nrow(mesh_dd))


for (i in 1:nrow(mesh_dd)){
  coll_vec = c() 
  for (j in (((mesh_dd$idx[i]-1)*num_sim)+1):(mesh_dd$idx[i]*num_sim)){
    tre = read.tree(paste0(tree_fp,"sim.",j,".tre"))
    char_dat =  read.table(paste0(tree_fp,'sim.',j,".extant.dat.csv"),sep = ",",header = T)
    # ignore tree with only 1 extant taxa 
    if (nrow(char_dat) != 1){
      tre_extant = drop.extinct(tre)
      coll_vec = append(coll_vec,colless_corr(tre_extant))
    }
  }
  mesh_dd$colless_norm[i] =mean(coll_vec,na.rm = TRUE)
}


# p1_ext = ggplot(mesh_ed, aes(x = e_d, y = colless_norm)) +
#   geom_line(color = "skyblue", linewidth = 1) +
#   labs(title = "normalized colless", x = expression(e^D), y = "normalized colless") +
#   scale_y_continuous(limits = c(0, 0.5))+
#   theme_minimal() + 
#   scale_x_reverse() 
# 
# p1_within = ggplot(mesh_wd, aes(x = w_d, y = colless_norm)) +
#   geom_line(color = "skyblue", linewidth = 1) +
#   labs(title = "normalized colless", x = expression(w^D), y = "normalized colless") +
#   scale_y_continuous(limits = c(0, 0.5))+
#   theme_minimal() 
# 
# p1_dispersal = ggplot(mesh_dd, aes(x = d_d_away, y = colless_norm)) +
#   geom_line(color = "skyblue", linewidth = 1) +
#   labs(title = "normalized colless", x = expression(d^{D*","*dest}), y = "normalized colless") +
#   scale_y_continuous(limits = c(0, 0.5))+
#   theme_minimal() 

p1_ext = ggplot(mesh_ed, aes(x = e_d, y = colless_norm)) +
  geom_line(color = "skyblue", linewidth = 1) +
  labs(title = "normalized colless", x = expression(e^D), y = "normalized colless") +
  scale_y_continuous(
    limits = c(0.008, 0.019),
    breaks = seq(0.008, 0.018, by = 0.002)
  ) +
  theme_minimal() + 
  scale_x_reverse() 

p1_within = ggplot(mesh_wd, aes(x = w_d, y = colless_norm)) +
  geom_line(color = "skyblue", linewidth = 1) +
  labs(title = "normalized colless", x = expression(w^D), y = "normalized colless") +
  scale_y_continuous(
    limits = c(0.008, 0.019),
    breaks = seq(0.008, 0.018, by = 0.002)
  )+
  theme_minimal() 

p1_dispersal = ggplot(mesh_dd, aes(x = d_d_away, y = colless_norm)) +
  geom_line(color = "skyblue", linewidth = 1) +
  labs(title = "normalized colless", x = expression(d^{D*","*dest}), y = "normalized colless") +
  scale_y_continuous(
    limits = c(0.008, 0.019),
    breaks = seq(0.008, 0.018, by = 0.002)
  )+
  theme_minimal() 

grid.arrange(arrangeGrob(p1_ext, p1_within, p1_dispersal,  p1_dispersal, nrow = 2, ncol = 2),heights = c(0.1, 1))

###################
# # e_d vs d_d_away heatmap (including w_d = 0 & != 0)
###################
p2_colless = ggplot(mesh, aes(x=d_d_away ,y=e_d,fill=colless_norm))
p2_colless = p2_colless + geom_tile() + 
  labs(title = expression("Effect of"~d^{D*","*dest}~"and"~e^D~"on"~"Colless"), x = expression(d^{D*","*dest}), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

p2_gamma = ggplot(mesh, aes(x=d_d_away ,y=e_d,fill=gamma))
p2_gamma = p2_gamma + geom_tile() + 
  labs(title = expression("Effect of"~d^{D*","*dest}~"and"~e^D~"on"~gamma), x = expression(d^{D*","*dest}), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

p2_size  = ggplot(mesh, aes(x=d_d_away ,y=e_d,fill=mean_range_size))
p2_size  = p2_size  + geom_tile() + 
  labs(title = expression("Effect of"~d^{D*","*dest}~"and"~e^D~"on mean range size"), x = expression(d^{D*","*dest}), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

p2_taxa  = ggplot(mesh, aes(x=d_d_away ,y=e_d,fill=log_mean_extant_taxa))
p2_taxa  = p2_taxa  + geom_tile() + 
  labs(title = expression("Effect of"~d^{D*","*dest}~"and"~e^D~"on log mean extant taxa"), x = expression(d^{D*","*dest}), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

###################
# # w_d vs e_d heatmap (including d_d_away = 0 & != 0)
###################
p1_colless = ggplot(mesh, aes(x=w_d ,y=e_d,fill=colless_norm))
p1_colless = p1_colless + geom_tile() + 
  labs(title = expression("Effect of"~w^D~"and"~e^D~"on"~"Colless"), x = expression(w^D), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

p1_gamma = ggplot(mesh, aes(x=w_d ,y=e_d,fill=gamma))
p1_gamma = p1_gamma + geom_tile() + 
  labs(title = expression("Effect of"~w^D~"and"~e^D~"on"~gamma), x = expression(w^D), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

p1_size  = ggplot(mesh, aes(x=w_d ,y=e_d,fill=mean_range_size))
p1_size  = p1_size  + geom_tile() + 
  labs(title = expression("Effect of"~w^D~"and"~e^D~"on mean range size"), x = expression(w^D), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

p1_taxa  = ggplot(mesh, aes(x=w_d ,y=e_d,fill=log_mean_extant_taxa))
p1_taxa  = p1_taxa  + geom_tile() + 
  labs(title = expression("Effect of"~w^D~"and"~e^D~"on log mean extant taxa"), x = expression(w^D), y = expression(e^D)) +
  scale_fill_viridis_c() +
  theme_minimal()

#### SAVE
# pdf(file = paste0(out_fp,"treestats_dd_colless.pdf"), width = 8, height = 6)
pdf(file = paste0(out_fp,"treestats_dd_small_colless.pdf"), width = 8, height = 6)

main_title <- textGrob(
  "Effect of diversity-dependence on",
  gp = gpar(fontsize = 16, fontface = "bold")
)
grid.arrange(main_title,arrangeGrob(p1_ext, p1_within, p1_dispersal,  p1_dispersal, nrow = 2, ncol = 2),heights = c(0.1, 1))


dev.off()


pdf(file = paste0(out_fp,"heatmap_wd_ed_colless.pdf"), width = 8, height = 6)
print(p1_colless)
dev.off()

pdf(file = paste0(out_fp,"heatmap_d_ed_colless.pdf"), width = 8, height = 6)
print(p2_colless)
dev.off()

pdf(file = paste0(out_fp,"heatmap_d_ed_gamma.pdf"), width = 8, height = 6)
print(p2_gamma)
dev.off()
