# This script compares sampled rate parameters with simulating distribution 
# (in case of base rates drawn from loguniform distribution)

fp      = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
in_fp   = paste0(fp,"data/deep_learning/server/DD-GeoSSE/5_reg/CURRENT_SET/model_selection_7/")
out_fp  = paste0(fp,"plot/") #output 

# model_types = c("ddall_7","dddisp_7","ddext_7","ddwithin_7",
#                 "ddext_disp_7","ddext_within_7","ddwithin_disp_7","ddgeosse_7"
#                 )
model_types = c("ddext_disp_7","ddext_within_7","ddwithin_disp_7","ddgeosse_7")

num_trees = 25000
# num_trees = 100

rho_w_true = exp(runif(num_trees, log(0.0001), log(1)))
rho_e_true = exp(runif(num_trees, log(0.0001), log(1)))
rho_d_true = exp(runif(num_trees, log(0.0001), log(1)))
rho_b_true = exp(runif(num_trees, log(0.0001), log(1)))
w_d_true   = runif(num_trees, -2, -0.0001)
d_d_true   = runif(num_trees, -2, -0.0001)
e_d_true   = runif(num_trees, 0.0001, 2)

for (i in 1:length(model_types)){
  par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
  rho_w_vec = c()
  rho_e_vec = c()
  rho_d_vec = c()
  rho_b_vec = c()
  #
  w_d_vec   = c()
  e_d_vec   = c()
  d_d_vec   = c()
  for (j in 1:num_trees){
    dat_par_j = read.csv(paste0(in_fp,model_types[i],"/simulate/sim.",j,".labels.csv"))
    #
    rho_w = exp(dat_par_j$log_base_w)
    rho_e = exp(dat_par_j$log_base_e)
    rho_d = exp(dat_par_j$log_base_d)
    rho_b = exp(dat_par_j$log_base_b)
    #
    w_d   = log(dat_par_j$exp_scalar_w)
    e_d   = log(dat_par_j$exp_scalar_e)
    d_d   = log(dat_par_j$exp_scalar_d_j)
    #
    rho_w_vec = append(rho_w_vec,rho_w)
    rho_e_vec = append(rho_e_vec,rho_e)
    rho_d_vec = append(rho_d_vec,rho_d)
    rho_b_vec = append(rho_b_vec,rho_b)
    #
    w_d_vec   = append(w_d_vec,w_d)
    e_d_vec   = append(e_d_vec,e_d)
    d_d_vec   = append(d_d_vec,d_d)
  }
  plot(density(rho_w_true), main = paste("rho_w:", model_types[i]), lwd = 2)
  lines(density(rho_w_vec), col = "blue", lwd = 2)
  
  plot(density(rho_e_true), main = paste("rho_e:", model_types[i]), lwd = 2)
  lines(density(rho_e_vec), col = "blue", lwd = 2)
  
  plot(density(rho_d_true), main = paste("rho_d:", model_types[i]), lwd = 2)
  lines(density(rho_d_vec), col = "blue", lwd = 2)
  
  plot(density(rho_b_true), main = paste("rho_b:", model_types[i]), lwd = 2)
  lines(density(rho_b_vec), col = "blue", lwd = 2)
  #
  if (model_types[i]=="ddall_7"){
    plot(density(w_d_true), main = paste("w_d:", model_types[i]), lwd = 2)
    lines(density(w_d_vec), col = "blue", lwd = 2)
    
    plot(density(e_d_true), main = paste("e_d:", model_types[i]), lwd = 2)
    lines(density(e_d_vec), col = "blue", lwd = 2)
    
    plot(density(d_d_true), main = paste("d_d:", model_types[i]), lwd = 2)
    lines(density(d_d_vec), col = "blue", lwd = 2)
  } else if (model_types[i]=="dddisp_7"){
    plot(density(d_d_true), main = paste("d_d:", model_types[i]), lwd = 2)
    lines(density(d_d_vec), col = "blue", lwd = 2)
  } else if (model_types[i]=="ddext_7"){
    plot(density(e_d_true), main = paste("e_d:", model_types[i]), lwd = 2)
    lines(density(e_d_vec), col = "blue", lwd = 2)
  } else if (model_types[i]=="ddwithin_7"){
    plot(density(w_d_true), main = paste("w_d:", model_types[i]), lwd = 2)
    lines(density(w_d_vec), col = "blue", lwd = 2)
  } else if (model_types[i]=="ddext_disp_7" ){
    plot(density(e_d_true), main = paste("e_d:", model_types[i]), lwd = 2)
    lines(density(e_d_vec), col = "blue", lwd = 2)
    
    plot(density(d_d_true), main = paste("d_d:", model_types[i]), lwd = 2)
    lines(density(d_d_vec), col = "blue", lwd = 2)
  } else if (model_types[i]=="ddext_within_7" ){
    plot(density(e_d_true), main = paste("e_d:", model_types[i]), lwd = 2)
    lines(density(e_d_vec), col = "blue", lwd = 2)
    
    plot(density(w_d_true), main = paste("w_d:", model_types[i]), lwd = 2)
    lines(density(w_d_vec), col = "blue", lwd = 2)
  } else if (model_types[i]=="ddwithin_disp_7" ){
    plot(density(w_d_true), main = paste("w_d:", model_types[i]), lwd = 2)
    lines(density(w_d_vec), col = "blue", lwd = 2)
    
    plot(density(d_d_true), main = paste("d_d:", model_types[i]), lwd = 2)
    lines(density(d_d_vec), col = "blue", lwd = 2)
  }
  #
  plot.new()
}