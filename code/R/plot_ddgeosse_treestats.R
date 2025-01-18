library(viridis)
library(plotly)
library(ggplot2)
library(gridExtra)
library(grid)

# Note: case where w_d = d_d_away = e_d = 0 (i.e. standard geosse) is excluded
###################
# # SET FILE SYSTEM
###################
fp              = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
out_fp          = paste0(fp,"plot/sim_treestats/")

###################
# # READ FILE
###################
mesh_dat = read.table(paste0(fp,"data/tree_stats/mesh_tree_stats_50.csv"),sep = ",",header = T) #based on extant-only trees
mesh_dat$log_mean_extant_taxa <- log(mesh_dat$mean_extant_taxa)

# Only dd-extinction
ext_dat = subset(mesh_dat,w_d == 0 & d_d_away == 0)

# (I think, better also to create barplots for each of these individual cases where density-dep is only applied to one event)

# Only dd-within 
within_dat = subset(mesh_dat,e_d == 0 & d_d_away == 0)
within_dat = within_dat[order(within_dat$w_d,decreasing = T), ]

# Only dd-dispersal
disperse_dat = subset(mesh_dat,w_d == 0 & e_d == 0)
disperse_dat = disperse_dat[order(disperse_dat$d_d_away,decreasing = T), ]

# Only dd-extinction & dispersal 
extdisperse_dat = subset(mesh_dat,w_d == 0 & d_d_away != 0 & e_d !=0)

# Only dd-within & dispersal
withindisperse_dat = subset(mesh_dat,e_d == 0 & w_d != 0 & d_d_away != 0)
withindisperse_dat = withindisperse_dat[order(withindisperse_dat$w_d,decreasing = T), ]

###################
# # HEATMAP PLOTTING
###################

###################
# # w_d vs e_d heatmap (including d_d_away = 0 & != 0)
###################
p1_beta = ggplot(mesh_dat, aes(x=w_d ,y=e_d,fill=beta))
p1_beta = p1_beta + geom_tile() + 
          labs(title = expression("Effect of"~w[d]~"and"~e[d]~"on"~beta), x = expression(w[d]), y = expression(e[d])) +
          scale_fill_viridis_c() +
          theme_minimal()

p1_gamma = ggplot(mesh_dat, aes(x=w_d ,y=e_d,fill=gamma))
p1_gamma = p1_gamma + geom_tile() + 
           labs(title = expression("Effect of"~w[d]~"and"~e[d]~"on"~gamma), x = expression(w[d]), y = expression(e[d])) +
           scale_fill_viridis_c() +
           theme_minimal()

p1_size  = ggplot(mesh_dat, aes(x=w_d ,y=e_d,fill=mean_range_size))
p1_size  = p1_size  + geom_tile() + 
           labs(title = expression("Effect of"~w[d]~"and"~e[d]~"on mean range size"), x = expression(w[d]), y = expression(e[d])) +
           scale_fill_viridis_c() +
           theme_minimal()

p1_taxa  = ggplot(mesh_dat, aes(x=w_d ,y=e_d,fill=log_mean_extant_taxa))
p1_taxa  = p1_taxa  + geom_tile() + 
           labs(title = expression("Effect of"~w[d]~"and"~e[d]~"on log mean extant taxa"), x = expression(w[d]), y = expression(e[d])) +
           scale_fill_viridis_c() +
           theme_minimal()

###################
# # e_d vs d_d_away heatmap (including w_d = 0 & != 0)
###################
p2_beta = ggplot(mesh_dat, aes(x=d_d_away ,y=e_d,fill=beta))
p2_beta = p2_beta + geom_tile() + 
          labs(title = expression("Effect of"~d[D*","*dest]~"and"~e[d]~"on"~beta), x = expression(d[D*","*dest]), y = expression(e[d])) +
          scale_fill_viridis_c() +
          theme_minimal()

p2_gamma = ggplot(mesh_dat, aes(x=d_d_away ,y=e_d,fill=gamma))
p2_gamma = p2_gamma + geom_tile() + 
           labs(title = expression("Effect of"~d[D*","*dest]~"and"~e[d]~"on"~gamma), x = expression(d[D*","*dest]), y = expression(e[d])) +
           scale_fill_viridis_c() +
           theme_minimal()

p2_size  = ggplot(mesh_dat, aes(x=d_d_away ,y=e_d,fill=mean_range_size))
p2_size  = p2_size  + geom_tile() + 
           labs(title = expression("Effect of"~d[D*","*dest]~"and"~e[d]~"on mean range size"), x = expression(d[D*","*dest]), y = expression(e[d])) +
           scale_fill_viridis_c() +
           theme_minimal()

p2_taxa  = ggplot(mesh_dat, aes(x=d_d_away ,y=e_d,fill=log_mean_extant_taxa))
p2_taxa  = p2_taxa  + geom_tile() + 
           labs(title = expression("Effect of"~d[D*","*dest]~"and"~e[d]~"on log mean extant taxa"), x = expression(d[D*","*dest]), y = expression(e[d])) +
           scale_fill_viridis_c() +
           theme_minimal()

###################
# # w_d vs d_d_away heatmap (including e_d = 0 & != 0)
###################
p3_beta = ggplot(mesh_dat, aes(x=w_d ,y=d_d_away,fill=beta))
p3_beta = p3_beta + geom_tile() + 
          labs(title = expression("Effect of"~d[D*","*dest]~"and"~w[d]~"on"~beta), x = expression(w[d]), y = expression(d[D*","*dest])) +
          scale_fill_viridis_c() +
          theme_minimal()

p3_gamma = ggplot(mesh_dat, aes(x=w_d ,y=d_d_away,fill=gamma))
p3_gamma = p3_gamma + geom_tile() + 
           labs(title = expression("Effect of"~d[D*","*dest]~"and"~w[d]~"on"~gamma), x = expression(w[d]), y = expression(d[D*","*dest])) +
           scale_fill_viridis_c() +
           theme_minimal()

p3_size  = ggplot(mesh_dat, aes(x=w_d ,y=d_d_away,fill=mean_range_size))
p3_size  = p3_size  + geom_tile() + 
           labs(title = expression("Effect of"~d[D*","*dest]~"and"~w[d]~"on mean range size"), x = expression(w[d]), y = expression(d[D*","*dest])) +
           scale_fill_viridis_c() +
           theme_minimal()

p3_taxa  = ggplot(mesh_dat, aes(x=w_d ,y=d_d_away,fill=log_mean_extant_taxa))
p3_taxa  = p3_taxa  + geom_tile() + 
           labs(title = expression("Effect of"~d[D*","*dest]~"and"~w[d]~"on log mean extant taxa"), x = expression(w[d]), y = expression(d[D*","*dest])) +
           scale_fill_viridis_c() +
           theme_minimal()

###################
# # BOXPLOT PLOTTING
###################

###################
# # density extinction
###################
p1_ext = ggplot(ext_dat, aes(x = e_d, y = beta)) +
  geom_line(color = "skyblue", size = 1) +
  labs(title = expression(beta), x = expression(e[d]), y = expression(beta)) +
  theme_minimal()

p2_ext = ggplot(ext_dat, aes(x = e_d, y = gamma)) +
  geom_line(color = "lightgreen",size = 1) +
  labs(title = expression(gamma), x = expression(e[d]), y = expression(gamma)) +
  theme_minimal()

p3_ext = ggplot(ext_dat, aes(x = e_d, y = mean_range_size)) +
  geom_line(color = "red", size = 1) +
  labs(title = expression("mean range size"), x = expression(e[d]), y = "mean range size") +
  theme_minimal()

p4_ext = ggplot(ext_dat, aes(x = e_d, y = mean_extant_taxa)) +
  geom_line(color = "black", size = 1) +
  labs(title = expression("mean extant taxa"), x = expression(e[d]), y = "mean extant taxa") +
  theme_minimal()

main_title_ext <- textGrob(
  "Effect of density-dependent extinction on",
  gp = gpar(fontsize = 16, fontface = "bold")
)

###################
# # density within
###################
p1_within = ggplot(within_dat, aes(x = w_d, y = beta)) +
  geom_line(color = "skyblue", size = 1) +
  labs(title = expression( beta), x = expression(w[d]), y = expression(beta)) +
  theme_minimal()

p2_within = ggplot(within_dat, aes(x = w_d, y = gamma)) +
  geom_line(color = "lightgreen", size = 1) +
  labs(title = expression(gamma), x = expression(w[d]), y = expression(gamma)) +
  theme_minimal()

p3_within = ggplot(within_dat, aes(x = w_d, y = mean_range_size)) +
  geom_line(color = "red", size = 1) +
  labs(title = expression("mean range size"), x = expression(w[d]), y = "mean range size") +
  theme_minimal()

p4_within = ggplot(within_dat, aes(x = w_d, y = mean_extant_taxa)) +
  geom_line(color = "black", size = 1) +
  labs(title = expression("mean extant taxa"), x = expression(w[d]), y = "mean extant taxa") +
  theme_minimal()

main_title_within <- textGrob(
  "Effect of density-dependent within-region speciation on",
  gp = gpar(fontsize = 16, fontface = "bold")
)

###################
# # density dispersal
###################
p1_dispersal = ggplot(disperse_dat, aes(x = d_d_away, y = beta)) +
  geom_line(color = "skyblue", size = 1) +
  labs(title = expression(beta), x = expression(d[D*","*dest]), y = expression(beta)) +
  theme_minimal()

p2_dispersal = ggplot(disperse_dat, aes(x = d_d_away, y = gamma)) +
  geom_line(color = "lightgreen", size = 1) +
  labs(title = expression(gamma), x = expression(d[D*","*dest]), y = expression(gamma)) +
  theme_minimal()

p3_dispersal = ggplot(disperse_dat, aes(x = d_d_away, y = mean_range_size)) +
  geom_line(color = "red", size = 1) +
  labs(title = expression("mean range size"), x = expression(d[D*","*dest]), y = "mean range size") +
  theme_minimal()

p4_dispersal = ggplot(disperse_dat, aes(x = d_d_away, y = mean_extant_taxa)) +
  geom_line(color = "black", size = 1) +
  labs(title = expression("mean extant taxa"), x = expression(d[D*","*dest]), y = "mean extant taxa") +
  theme_minimal()

main_title_disperse <- textGrob(
  "Effect of density-dependent incoming dispersal on",
  gp = gpar(fontsize = 16, fontface = "bold")
)

###################
# # SAVE
###################
pdf(file = paste0(out_fp,"treestats_dd.pdf"), width = 8, height = 6)

# Print each plot to a new page in the PDF
print(p1_beta)
print(p1_gamma)
print(p1_size)
print(p1_taxa)
#
print(p2_beta)
print(p2_gamma)
print(p2_size)
print(p2_taxa)
#
print(p3_beta)
print(p3_gamma)
print(p3_size)
print(p3_taxa)
#
grid.arrange(main_title_ext,arrangeGrob(p1_ext, p2_ext, p3_ext, p4_ext, nrow = 2, ncol = 2),heights = c(0.1, 1))
grid.arrange(main_title_within,arrangeGrob(p1_within, p2_within, p3_within, p4_within ,nrow = 2,ncol = 2),heights = c(0.1, 1))
grid.arrange(main_title_disperse,arrangeGrob(p1_dispersal, p2_dispersal, p3_dispersal, p4_dispersal,nrow = 2, ncol = 2),heights = c(0.1, 1) )

#
dev.off()
