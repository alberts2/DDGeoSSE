# This figure plots returns carrying capacity for parameters from Figure 4
library(nleqslv)

fp      = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
# tree_fp = paste0(fp,"data/tree_stats/smaller_interval/data_treestats_small/")
tree_fp = paste0(fp,"data/tree_stats/larger_interval/data_treestats_big/")
# dat_fp  = paste0(fp,"data/tree_stats/smaller_interval/")
dat_fp  = paste0(fp,"data/tree_stats/larger_interval/")
out_fp          = paste0(fp,"plot/sim_treestats/")

num_sim = 50 
# mesh = read.table(paste0(dat_fp,'mesh_tree_stats_50_small.csv'),sep = ",",header = T)
mesh = read.table(paste0(dat_fp,'mesh_tree_stats_50_big.csv'),sep = ",",header = T)
mesh$log_n_A = rep(NA,nrow(mesh))
mesh$log_n_B = rep(NA,nrow(mesh))
mesh$log_n_C = rep(NA,nrow(mesh))

for (i in 1:nrow(mesh)){
  # Give set of parameters, find n_i_start for each region under unconstrained model
  
  rho_w = 0.05 #base rate for speciation
  rho_d = 0.05 #base rate for dispersal
  rho_e = 0.02 #base rate for extinction
  
  
  wD      = mesh$w_d[i] #DD effect for speciation
  dD_dest = mesh$d_d_away[i] #DD effect for inbound dispersal
  dD_src  = 0 #DD effect for outbound dispersal
  eD      = mesh$e_d[i] #DD effect for extinction
  
  
  # use log_n for the search so it does not fail if it searches negative values
  equil_fun_log <- function(log_n) {
    n   = exp(log_n) #vector of ni* to find
    L   = log(n) + 1
    L_2 = log(n+1)+1
    K   = length(n)
    
    F = numeric(K)
    
    # balance equation for each region
    for (i in seq_len(K)) {
      incoming = sum(L[-i]^dD_src)
      
      F[i] = rho_w * L[i]^wD + rho_d * L_2[i]^dD_dest * incoming -rho_e * L[i]^eD
    }
    F
  }
  #
  fit <- nleqslv(
    x = log(c(1, 1, 1)), #starting value for each ni_star
    fn = equil_fun_log
  )
  mesh$log_n_A[i] = format(round(fit$x[1], 0), scientific = FALSE)
  mesh$log_n_B[i] = format(round(fit$x[2], 0), scientific = FALSE)
  mesh$log_n_C[i] = format(round(fit$x[3], 0), scientific = FALSE)
}

ext=mesh[mesh$d_d_away==0 & mesh$e_d!=0 & mesh$w_d==0,]
within=mesh[mesh$d_d_away==0 & mesh$e_d==0 & mesh$w_d!=0,]
dispersal=mesh[mesh$d_d_away!=0 & mesh$e_d==0 & mesh$w_d==0,]

p1_ext = ggplot(ext, aes(x = e_d, y = as.numeric(log_n_A))) +
  geom_line(color = "skyblue", linewidth = 1) +
  labs(title = "ddext", x = expression(e^D), y = "log n_i") +
  theme_minimal() + 
  scale_y_log10()+
  scale_x_reverse() 

p1_within = ggplot(within, aes(x = w_d, y = as.numeric(log_n_A))) +
  geom_line(color = "skyblue", linewidth = 1) +
  labs(title = "ddwithin", x = expression(w^D), y = "log n_i") +
  theme_minimal() + 
  scale_y_log10()+
  scale_x_reverse() 