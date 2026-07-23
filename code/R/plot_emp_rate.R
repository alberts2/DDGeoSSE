# This script plots process rates using estimated parameters then
# overlay point estimate of empirical rate for each region on the curve 
library(ape)
library(ggplot2)

###### FILE SETTING ######
fp              = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
anolis_fp       = paste0(fp,"data/empiric_anolis/use_this/")
viburnum_fp     = paste0(fp,"data/empiric_viburnum_rescaled/")
# empirical estimate for anolis from trained network 1 (no_weight)
est_anolis_1_fp = paste0(fp,"data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/parameter_estimation/Anolis/final_submodel_sum_no_weight/submodel_7_mixture_anolis/estimate/")
# empirical estimate for viburnum from trained network 1 (no_weight)
est_viburnum_1_fp = paste0(fp,"data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/parameter_estimation/Viburnum/final_submodel_sum_no_weight/submodel_6_mixture_oreino/estimate/")
# empirical estimate for anolis from trained network 2 (weight)
est_anolis_2_fp = paste0(fp,"data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/parameter_estimation/Anolis/final_submodel_sum_weight_0.001/submodel_7_mixture_anolis/estimate/")
# empirical estimate for viburnum from trained network 2 (weight)
est_viburnum_2_fp = paste0(fp,"data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/parameter_estimation/Viburnum/final_submodel_sum_weight_0.001/submodel_6_mixture_oreino/estimate/")
#
out_fp          = paste0(fp,"plot/emp_rate/")

setwd(fp)

# load anolis tree
anolis_tre = read.tree(paste0(anolis_fp,"anolis.5.tre"))
# load viburnum tree
viburnum_tre = read.tree(paste0(viburnum_fp,"oreino5.1.tre"))
# load anolis' character data
anolis_char = read.table(paste0(anolis_fp,"anolis.5.dat.csv"),sep = ",",header = T)
# load viburnum' character data
viburnum_char = read.table(paste0(viburnum_fp,"oreino5.1.dat.csv"),sep = ",",header = T)

# load empirical estimate for anolis (network 1)
est_anolis_1 = read.table(paste0(est_anolis_1_fp,"sim.empirical_est.labels_num.csv"),sep = ",",header = T)
anolis_base_vec_1 = c(exp(est_anolis_1$log_base_w_value),
                      exp(est_anolis_1$log_base_e_value),
                      exp(est_anolis_1$log_base_d_value)) # for w,e and d (linear scale)
#
anolis_lower_vec_1 = c(exp(est_anolis_1$log_base_w_lower),
                       exp(est_anolis_1$log_base_e_lower),
                       exp(est_anolis_1$log_base_d_lower)) # for lower w,e and d (linear scale)
#
anolis_upper_vec_1 = c(exp(est_anolis_1$log_base_w_upper),
                       exp(est_anolis_1$log_base_e_upper),
                       exp(est_anolis_1$log_base_d_upper)) # for upper w,e and d (linear scale)
#
anolis_dd_vec_1  = c(log(est_anolis_1$exp_scalar_w_value),
                     log(est_anolis_1$exp_scalar_e_value),
                     log(est_anolis_1$exp_scalar_d_j_value)) #DD scalar for w, e and incoming d (linear scale)
#
anolis_lower_dd_vec_1  = c(log(est_anolis_1$exp_scalar_w_lower),
                           log(est_anolis_1$exp_scalar_e_lower),
                           log(est_anolis_1$exp_scalar_d_j_lower)) # lower DD scalar for w, e and incoming d (linear scale)
#
anolis_upper_dd_vec_1  = c(log(est_anolis_1$exp_scalar_w_upper),
                           log(est_anolis_1$exp_scalar_e_upper),
                           log(est_anolis_1$exp_scalar_d_j_upper)) # upper DD scalar for w, e and incoming d (linear scale)

# load empirical estimate for anolis (network 2)
est_anolis_2 = read.table(paste0(est_anolis_2_fp,"sim.empirical_est.labels_num.csv"),sep = ",",header = T)
anolis_base_vec_2 = c(exp(est_anolis_2$log_base_w_value),
                      exp(est_anolis_2$log_base_e_value),
                      exp(est_anolis_2$log_base_d_value)) # for w,e and d (linear scale)
#
anolis_lower_vec_2 = c(exp(est_anolis_2$log_base_w_lower),
                       exp(est_anolis_2$log_base_e_lower),
                       exp(est_anolis_2$log_base_d_lower)) # for lower w,e and d (linear scale)
#
anolis_upper_vec_2 = c(exp(est_anolis_2$log_base_w_upper),
                       exp(est_anolis_2$log_base_e_upper),
                       exp(est_anolis_2$log_base_d_upper)) # for upper w,e and d (linear scale)
#
anolis_dd_vec_2  = c(log(est_anolis_2$exp_scalar_w_value),
                     log(est_anolis_2$exp_scalar_e_value),
                     log(est_anolis_2$exp_scalar_d_j_value)) #DD scalar for w, e and incoming d (linear scale)
#
anolis_lower_dd_vec_2  = c(log(est_anolis_2$exp_scalar_w_lower),
                           log(est_anolis_2$exp_scalar_e_lower),
                           log(est_anolis_2$exp_scalar_d_j_lower)) # lower DD scalar for w, e and incoming d (linear scale)
#
anolis_upper_dd_vec_2  = c(log(est_anolis_2$exp_scalar_w_upper),
                           log(est_anolis_2$exp_scalar_e_upper),
                           log(est_anolis_2$exp_scalar_d_j_upper)) # upper DD scalar for w, e and incoming d (linear scale)
#
# load empirical estimate for viburnum (network 1)
est_viburnum_1 = read.table(paste0(est_viburnum_1_fp,"sim.empirical_est.labels_num.csv"),sep = ",",header = T)
viburnum_base_vec_1 = c(7*exp(est_viburnum_1$log_base_w_value),
                        7*exp(est_viburnum_1$log_base_d_value)) # for w and d (linear scale)
#
viburnum_lower_vec_1 = c(7*exp(est_viburnum_1$log_base_w_lower),
                         7*exp(est_viburnum_1$log_base_d_lower)) # for lower w and d (linear scale)
#
viburnum_upper_vec_1 = c(7*exp(est_viburnum_1$log_base_w_upper),
                         7*exp(est_viburnum_1$log_base_d_upper)) # for upper w andd (linear scale)
#
viburnum_dd_vec_1  = c(log(est_viburnum_1$exp_scalar_w_value),
                       log(est_viburnum_1$exp_scalar_d_j_value)) #DD scalar for w and incoming d (linear scale)
#
viburnum_lower_dd_vec_1  = c(log(est_viburnum_1$exp_scalar_w_lower),
                             log(est_viburnum_1$exp_scalar_d_j_lower)) # lower DD scalar for w and incoming d (linear scale)
#
viburnum_upper_dd_vec_1  = c(log(est_viburnum_1$exp_scalar_w_upper),
                             log(est_viburnum_1$exp_scalar_d_j_upper)) # upper DD scalar for w and incoming d (linear scale)
#
# load empirical estimate for viburnum (network 2)
est_viburnum_2 = read.table(paste0(est_viburnum_2_fp,"sim.empirical_est.labels_num.csv"),sep = ",",header = T)
viburnum_base_vec_2 = c(7*exp(est_viburnum_2$log_base_w_value),
                        7*exp(est_viburnum_2$log_base_d_value)) # for w and d (linear scale)
#
viburnum_lower_vec_2 = c(7*exp(est_viburnum_2$log_base_w_lower),
                         7*exp(est_viburnum_2$log_base_d_lower)) # for lower w and d (linear scale)
#
viburnum_upper_vec_2 = c(7*exp(est_viburnum_2$log_base_w_upper),
                         7*exp(est_viburnum_2$log_base_d_upper)) # for upper w and d (linear scale)
#
viburnum_dd_vec_2  = c(log(est_viburnum_2$exp_scalar_w_value),
                       log(est_viburnum_2$exp_scalar_d_j_value)) #DD scalar for w and incoming d (linear scale)
#
viburnum_lower_dd_vec_2  = c(log(est_viburnum_2$exp_scalar_w_lower),
                             log(est_viburnum_2$exp_scalar_d_j_lower)) # lower DD scalar for w and incoming d (linear scale)
#
viburnum_upper_dd_vec_2  = c(log(est_viburnum_2$exp_scalar_w_upper),
                             log(est_viburnum_2$exp_scalar_d_j_upper)) # upper DD scalar for w and incoming d (linear scale)


# count number of extant species in each region (anolis)
num_LA = sum(anolis_char[,2]) #Lesser Antilles
num_PR = sum(anolis_char[,3]) #Puerto Rico
num_CU = sum(anolis_char[,4]) #Cuba + Caymans
num_HI = sum(anolis_char[,5]) #Hispaniola
num_JA = sum(anolis_char[,6]) #Jamaica
num_spec_anolis_vec = c(num_LA,num_PR,num_CU,num_HI,num_JA)

# count number of extant species in each region (viburnum)
num_MX = sum(viburnum_char[,2]) #E/C/W Mexico
num_OA = sum(viburnum_char[,3]) #Oaxaca + Chiapas + Upper CA
num_LC = sum(viburnum_char[,4]) #Lower CA 
num_NA = sum(viburnum_char[,5]) #Northern Andes
num_CA = sum(viburnum_char[,6]) #Central Andes
num_spec_viburnum_vec = c(num_MX,num_OA,num_LC,num_NA,num_CA)

# define function that compute tip rates for all process 
compute_rate_anolis <- function(base_vec,dd_vec,num_spec,num_reg){
  # compute within-region speciation
  base_w = base_vec[1]
  w_D    = dd_vec[1]
  w_rate = base_w*(log(num_spec)+1)^w_D
  # w_rate = (log(num_spec)+1)^w_D
  # compute extinction
  base_e = base_vec[2]
  e_D    = dd_vec[2]
  e_rate = base_e*(log(num_spec)+1)^e_D
  # e_rate = (log(num_spec)+1)^e_D
  # compute total incoming dispersal
  # d_j = sum(d_ij) 
  # note d_D^src = 0
  base_d = base_vec[3]
  d_D_dest = dd_vec[3]
  d_rate = (num_reg-1)*base_d*1*(log(num_spec+1)+1)^d_D_dest
  # d_rate = (num_reg-1)*1*(log(num_spec+1)+1)^d_D_dest
  #
  return(list(
    w_rate,
    e_rate,
    d_rate
  ))
}

compute_rate_viburnum <- function(base_vec,dd_vec,num_spec,num_reg){
  # compute within-region speciation
  base_w = base_vec[1]
  w_D    = dd_vec[1]
  w_rate = base_w*(log(num_spec)+1)^w_D
  # w_rate = (log(num_spec)+1)^w_D
  # compute total incoming dispersal
  # d_j = sum(d_ij) 
  # note d_D^src = 0
  base_d = base_vec[2]
  d_D_dest = dd_vec[2]
  d_rate = (num_reg-1)*base_d*1*(log(num_spec+1)+1)^d_D_dest
  # d_rate = (num_reg-1)*1*(log(num_spec+1)+1)^d_D_dest
  #
  return(list(
    w_rate,
    d_rate
  ))
}

n_seq <- seq(1, 200, by = 1) # number of species in a region

# get rates for anolis from network 1 
curve_anolis_df_1 <- data.frame(
  num_spec = n_seq,
  within_region_speciation = rep(NA,200),
  within_region_speciation_low = rep(NA,200),
  within_region_speciation_high = rep(NA,200),
  extinction = rep(NA,200),
  extinction_low = rep(NA,200),
  extinction_high = rep(NA,200),
  total_in_disp = rep(NA,200),
  total_in_disp_low = rep(NA,200),
  total_in_disp_high = rep(NA,200)
)

# get rates for anolis from network 2 
curve_anolis_df_2 <- data.frame(
  num_spec = n_seq,
  within_region_speciation = rep(NA,200),
  within_region_speciation_low = rep(NA,200),
  within_region_speciation_high = rep(NA,200),
  extinction = rep(NA,200),
  extinction_low = rep(NA,200),
  extinction_high = rep(NA,200),
  total_in_disp = rep(NA,200),
  total_in_disp_low = rep(NA,200),
  total_in_disp_high = rep(NA,200)
)

# get rates for viburnum from network 1
curve_viburnum_df_1 <- data.frame(
  num_spec = n_seq,
  within_region_speciation = rep(NA,200),
  within_region_speciation_low = rep(NA,200),
  within_region_speciation_high = rep(NA,200),
  total_in_disp = rep(NA,200),
  total_in_disp_low = rep(NA,200),
  total_in_disp_high = rep(NA,200)
)

# get rates for viburnum from network 2
curve_viburnum_df_2 <- data.frame(
  num_spec = n_seq,
  within_region_speciation = rep(NA,200),
  within_region_speciation_low = rep(NA,200),
  within_region_speciation_high = rep(NA,200),
  total_in_disp = rep(NA,200),
  total_in_disp_low = rep(NA,200),
  total_in_disp_high = rep(NA,200)
)




for (i in 1:length(n_seq)){
  # compute the rate for anolis from network 1 for each n 
  get_rate_anolis_1 = compute_rate_anolis(base_vec = anolis_base_vec_1, dd_vec = anolis_dd_vec_1,
                                          num_spec = n_seq[i],num_reg = 5)
  get_lower_rate_anolis_1 = compute_rate_anolis(base_vec = anolis_lower_vec_1, dd_vec = anolis_lower_dd_vec_1,
                                                num_spec = n_seq[i],num_reg = 5)
  get_upper_rate_anolis_1 = compute_rate_anolis(base_vec = anolis_upper_vec_1, dd_vec = anolis_upper_dd_vec_1,
                                                num_spec = n_seq[i],num_reg = 5)
  # compute the rate for anolis from network 2 for each n 
  get_rate_anolis_2 = compute_rate_anolis(base_vec = anolis_base_vec_2, dd_vec = anolis_dd_vec_2,
                                          num_spec = n_seq[i],num_reg = 5)
  get_lower_rate_anolis_2 = compute_rate_anolis(base_vec = anolis_lower_vec_2, dd_vec = anolis_lower_dd_vec_2,
                                                num_spec = n_seq[i],num_reg = 5)
  get_upper_rate_anolis_2 = compute_rate_anolis(base_vec = anolis_upper_vec_2, dd_vec = anolis_upper_dd_vec_2,
                                                num_spec = n_seq[i],num_reg = 5)
  # compute the rate for viburnum from network 1 for each n 
  get_rate_viburnum_1 = compute_rate_viburnum(base_vec = viburnum_base_vec_1, dd_vec = viburnum_dd_vec_1,
                                              num_spec = n_seq[i],num_reg = 5)
  get_lower_rate_viburnum_1 = compute_rate_viburnum(base_vec = viburnum_lower_vec_1, dd_vec = viburnum_lower_dd_vec_1,
                                                    num_spec = n_seq[i],num_reg = 5)
  get_upper_rate_viburnum_1 = compute_rate_viburnum(base_vec = viburnum_upper_vec_1, dd_vec = viburnum_upper_dd_vec_1,
                                                    num_spec = n_seq[i],num_reg = 5)
  # compute the rate for viburnum from network 2 for each n 
  get_rate_viburnum_2 = compute_rate_viburnum(base_vec = viburnum_base_vec_2, dd_vec = viburnum_dd_vec_2,
                                              num_spec = n_seq[i],num_reg = 5)
  get_lower_rate_viburnum_2 = compute_rate_viburnum(base_vec = viburnum_lower_vec_2, dd_vec = viburnum_lower_dd_vec_2,
                                                    num_spec = n_seq[i],num_reg = 5)
  get_upper_rate_viburnum_2 = compute_rate_viburnum(base_vec = viburnum_upper_vec_2, dd_vec = viburnum_upper_dd_vec_2,
                                                    num_spec = n_seq[i],num_reg = 5)
  # append to the dataframes
  curve_anolis_df_1[i,2:ncol(curve_anolis_df_1)] = c(get_rate_anolis_1[1],get_lower_rate_anolis_1[1],get_upper_rate_anolis_1[1],
                                                     get_rate_anolis_1[2],get_lower_rate_anolis_1[2],get_upper_rate_anolis_1[2],
                                                     get_rate_anolis_1[3],get_lower_rate_anolis_1[3],get_upper_rate_anolis_1[3])
  #
  curve_anolis_df_2[i,2:ncol(curve_anolis_df_2)] = c(get_rate_anolis_2[1],get_lower_rate_anolis_2[1],get_upper_rate_anolis_2[1],
                                                     get_rate_anolis_2[2],get_lower_rate_anolis_2[2],get_upper_rate_anolis_2[2],
                                                     get_rate_anolis_2[3],get_lower_rate_anolis_2[3],get_upper_rate_anolis_2[3])
  #
  curve_viburnum_df_1[i,2:ncol(curve_viburnum_df_1)] = c(get_rate_viburnum_1[1],get_lower_rate_viburnum_1[1],get_upper_rate_viburnum_1[1],
                                                         get_rate_viburnum_1[2],get_lower_rate_viburnum_1[2],get_upper_rate_viburnum_1[2])
  #
  curve_viburnum_df_2[i,2:ncol(curve_viburnum_df_2)] = c(get_rate_viburnum_2[1],get_lower_rate_viburnum_2[1],get_upper_rate_viburnum_2[1],
                                                         get_rate_viburnum_2[2],get_lower_rate_viburnum_2[2],get_upper_rate_viburnum_2[2])
}

## get the actual rates for anolis and viburnum with their number of species in each region
## ANOLIS 
# anolis speciation rate for network 1 
anolis_speciation_dat_1 = curve_anolis_df_1[curve_anolis_df_1$num_spec %in% num_spec_anolis_vec,
                                            c("num_spec","within_region_speciation")]
# anolis speciation rate for network 2 
anolis_speciation_dat_2 = curve_anolis_df_2[curve_anolis_df_2$num_spec %in% num_spec_anolis_vec,
                                            c("num_spec","within_region_speciation")]
# anolis extinction rate for network 1
anolis_extinction_dat_1 = curve_anolis_df_1[curve_anolis_df_1$num_spec %in% num_spec_anolis_vec,
                                            c("num_spec","extinction")]
# anolis extinction rate for network 2
anolis_extinction_dat_2 = curve_anolis_df_2[curve_anolis_df_2$num_spec %in% num_spec_anolis_vec,
                                            c("num_spec","extinction")]
# anolis total dispersal rate for network 1
anolis_dispersal_dat_1 = curve_anolis_df_1[curve_anolis_df_1$num_spec %in% num_spec_anolis_vec,
                                           c("num_spec","total_in_disp")]
# anolis total dispersal rate for network 2
anolis_dispersal_dat_2 = curve_anolis_df_2[curve_anolis_df_2$num_spec %in% num_spec_anolis_vec,
                                           c("num_spec","total_in_disp")]

## Viburnum
# viburnum speciation rate for network 1 
idx <- match(num_spec_viburnum_vec, curve_viburnum_df_1$num_spec)

viburnum_speciation_dat_1 <- curve_viburnum_df_1[idx,c("num_spec", "within_region_speciation")]

# viburnum speciation rate for network 2 
idx_2 <- match(num_spec_viburnum_vec, curve_viburnum_df_2$num_spec)

viburnum_speciation_dat_2 <- curve_viburnum_df_2[idx_2,c("num_spec", "within_region_speciation")]

# viburnum total dispersal rate for network 1
viburnum_dispersal_dat_1 = curve_viburnum_df_1[idx,c("num_spec","total_in_disp")]
# viburnum total dispersal rate for network 2
viburnum_dispersal_dat_2 = curve_viburnum_df_2[idx_2,c("num_spec","total_in_disp")]

# Define regions and their level
## anolis
anolis_speciation_dat_1$region = c("Jamaica", "Puerto Rico", "Lesser Antilles", 
                                   "Hispaniola", "Cuba + Caymans")
anolis_speciation_dat_2$region = c("Jamaica", "Puerto Rico", "Lesser Antilles", 
                                   "Hispaniola", "Cuba + Caymans")
anolis_extinction_dat_1$region = c("Jamaica", "Puerto Rico", "Lesser Antilles", 
                                   "Hispaniola", "Cuba + Caymans")
anolis_extinction_dat_2$region = c("Jamaica", "Puerto Rico", "Lesser Antilles", 
                                   "Hispaniola", "Cuba + Caymans")
anolis_dispersal_dat_1$region = c("Jamaica", "Puerto Rico", "Lesser Antilles", 
                                  "Hispaniola", "Cuba + Caymans")
anolis_dispersal_dat_2$region = c("Jamaica", "Puerto Rico", "Lesser Antilles", 
                                  "Hispaniola", "Cuba + Caymans")

anolis_speciation_dat_1$region = factor(anolis_speciation_dat_1$region, levels = c("Lesser Antilles", "Puerto Rico", "Cuba + Caymans", 
                                                                                   "Hispaniola", "Jamaica"))
anolis_speciation_dat_2$region = factor(anolis_speciation_dat_2$region, levels = c("Lesser Antilles", "Puerto Rico", "Cuba + Caymans", 
                                                                                   "Hispaniola", "Jamaica"))
anolis_extinction_dat_1$region = factor(anolis_extinction_dat_1$region, levels = c("Lesser Antilles", "Puerto Rico", "Cuba + Caymans", 
                                                                                   "Hispaniola", "Jamaica"))
anolis_extinction_dat_2$region = factor(anolis_extinction_dat_2$region, levels = c("Lesser Antilles", "Puerto Rico", "Cuba + Caymans", 
                                                                                   "Hispaniola", "Jamaica"))
anolis_dispersal_dat_1$region = factor(anolis_dispersal_dat_1$region, levels = c("Lesser Antilles", "Puerto Rico", "Cuba + Caymans", 
                                                                                 "Hispaniola", "Jamaica"))
anolis_dispersal_dat_2$region = factor(anolis_dispersal_dat_2$region, levels = c("Lesser Antilles", "Puerto Rico", "Cuba + Caymans", 
                                                                                 "Hispaniola", "Jamaica"))
## viburnum
viburnum_speciation_dat_1$region = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                     "Lower CA","Northern Andes","Central Andes")
viburnum_speciation_dat_2$region = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                     "Lower CA","Northern Andes","Central Andes")
viburnum_dispersal_dat_1$region = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                    "Lower CA","Northern Andes","Central Andes")
viburnum_dispersal_dat_2$region = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                    "Lower CA","Northern Andes","Central Andes")

viburnum_speciation_dat_1$region = factor(viburnum_speciation_dat_1$region, levels = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                                                                       "Lower CA","Northern Andes","Central Andes"))
viburnum_speciation_dat_2$region = factor(viburnum_speciation_dat_2$region, levels = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                                                                       "Lower CA","Northern Andes","Central Andes"))
viburnum_dispersal_dat_1$region = factor(viburnum_dispersal_dat_1$region, levels = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                                                                     "Lower CA","Northern Andes","Central Andes"))
viburnum_dispersal_dat_2$region = factor(viburnum_dispersal_dat_2$region, levels = c("E/C/W Mexico","Oaxaca + Chiapas + Upper CA",
                                                                                     "Lower CA","Northern Andes","Central Andes"))

### PLOTTING
p_spec_anolis_1 = ggplot() +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = within_region_speciation_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = within_region_speciation_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = within_region_speciation),
            linewidth = 1) +
  geom_point(data = anolis_speciation_dat_1,
             aes(x = num_spec, y = within_region_speciation, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "Lesser Antilles" = "#d62728",  # red
    "Puerto Rico" = "#FFD700",      # yellow
    "Cuba + Caymans" = "#17becf",   # cyan
    "Hispaniola" = "#8c2bbd",       # purple
    "Jamaica" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "Lesser Antilles" = 1,   # hollow circle
    "Puerto Rico" = 17,  # triangle
    "Cuba + Caymans" = 15,       # square
    "Hispaniola" = 18, # diamond
    "Jamaica" = 3    # plus
  )) + 
  labs(
    x = "Num. extant species",
    y = "Within-region speciation rate"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_spec_anolis_2 = ggplot() +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = within_region_speciation_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = within_region_speciation_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = within_region_speciation),
            linewidth = 1) +
  geom_point(data = anolis_speciation_dat_2,
             aes(x = num_spec, y = within_region_speciation, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "Lesser Antilles" = "#d62728",  # red
    "Puerto Rico" = "#FFD700",      # yellow
    "Cuba + Caymans" = "#17becf",   # cyan
    "Hispaniola" = "#8c2bbd",       # purple
    "Jamaica" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "Lesser Antilles" = 1,   # hollow circle
    "Puerto Rico" = 17,  # triangle
    "Cuba + Caymans" = 15,       # square
    "Hispaniola" = 18, # diamond
    "Jamaica" = 3   # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Within-region speciation rate"
  ) +
  theme_classic() + 
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_ext_anolis_1 = ggplot() +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = extinction_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = extinction_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = extinction),
            linewidth = 1) +
  geom_point(data = anolis_extinction_dat_1,
             aes(x = num_spec, y = extinction, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "Lesser Antilles" = "#d62728",  # red
    "Puerto Rico" = "#FFD700",      # yellow
    "Cuba + Caymans" = "#17becf",   # cyan
    "Hispaniola" = "#8c2bbd",       # purple
    "Jamaica" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "Lesser Antilles" = 1,   # hollow circle
    "Puerto Rico" = 17,  # triangle
    "Cuba + Caymans" = 15,       # square
    "Hispaniola" = 18, # diamond
    "Jamaica" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Extinction rate"
  ) +
  theme_classic() + 
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_ext_anolis_2 = ggplot() +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = extinction_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = extinction_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = extinction),
            linewidth = 1) +
  geom_point(data = anolis_extinction_dat_2,
             aes(x = num_spec, y = extinction, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "Lesser Antilles" = "#d62728",  # red
    "Puerto Rico" = "#FFD700",      # yellow
    "Cuba + Caymans" = "#17becf",   # cyan
    "Hispaniola" = "#8c2bbd",       # purple
    "Jamaica" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "Lesser Antilles" = 1,   # hollow circle
    "Puerto Rico" = 17,  # triangle
    "Cuba + Caymans" = 15,       # square
    "Hispaniola" = 18, # diamond
    "Jamaica" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Extinction rate"
  ) +
  theme_classic() + 
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_disp_anolis_1 = ggplot() +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = total_in_disp_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = total_in_disp_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_1,
            aes(x = num_spec, y = total_in_disp),
            linewidth = 1) +
  geom_point(data = anolis_dispersal_dat_1,
             aes(x = num_spec, y = total_in_disp, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "Lesser Antilles" = "#d62728",  # red
    "Puerto Rico" = "#FFD700",      # yellow
    "Cuba + Caymans" = "#17becf",   # cyan
    "Hispaniola" = "#8c2bbd",       # purple
    "Jamaica" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "Lesser Antilles" = 1,   # hollow circle
    "Puerto Rico" = 17,  # triangle
    "Cuba + Caymans" = 15,       # square
    "Hispaniola" = 18, # diamond
    "Jamaica" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Total inbound dispersal rate"
  ) +
  theme_classic() + 
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_disp_anolis_2 = ggplot() +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = total_in_disp_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = total_in_disp_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_anolis_df_2,
            aes(x = num_spec, y = total_in_disp),
            linewidth = 1) +
  geom_point(data = anolis_dispersal_dat_2,
             aes(x = num_spec, y = total_in_disp, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "Lesser Antilles" = "#d62728",  # red
    "Puerto Rico" = "#FFD700",      # yellow
    "Cuba + Caymans" = "#17becf",   # cyan
    "Hispaniola" = "#8c2bbd",       # purple
    "Jamaica" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "Lesser Antilles" = 1,   # hollow circle
    "Puerto Rico" = 17,  # triangle
    "Cuba + Caymans" = 15,       # square
    "Hispaniola" = 18, # diamond
    "Jamaica" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Total inbound dispersal rate"
  ) +
  theme_classic() + 
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_spec_viburnum_1 = ggplot() +
  geom_line(data = curve_viburnum_df_1,
            aes(x = num_spec, y = within_region_speciation_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_1,
            aes(x = num_spec, y = within_region_speciation_low),
            linetype = "dashed",linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_1,
            aes(x = num_spec, y = within_region_speciation),
            linewidth = 1) +
  geom_point(data = viburnum_speciation_dat_1,
             aes(x = num_spec, y = within_region_speciation, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "E/C/W Mexico" = "#d62728",  # red
    "Oaxaca + Chiapas + Upper CA" = "#FFD700",      # yellow
    "Lower CA" = "#17becf",   # cyan
    "Northern Andes" = "#8c2bbd",       # purple
    "Central Andes" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "E/C/W Mexico" = 1,   # hollow circle
    "Oaxaca + Chiapas + Upper CA" = 17,  # triangle
    "Lower CA" = 15,       # square
    "Northern Andes" = 18, # diamond
    "Central Andes" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Within-region speciation rate"
  ) +
  theme_classic()+
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()


p_spec_viburnum_2 = ggplot() +
  geom_line(data = curve_viburnum_df_2,
            aes(x = num_spec, y = within_region_speciation_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_2,
            aes(x = num_spec, y = within_region_speciation_low),
            linetype = "dashed",linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_2,
            aes(x = num_spec, y = within_region_speciation),
            linewidth = 1) +
  geom_point(data = viburnum_speciation_dat_2,
             aes(x = num_spec, y = within_region_speciation, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "E/C/W Mexico" = "#d62728",  # red
    "Oaxaca + Chiapas + Upper CA" = "#FFD700",      # yellow
    "Lower CA" = "#17becf",   # cyan
    "Northern Andes" = "#8c2bbd",       # purple
    "Central Andes" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "E/C/W Mexico" = 1,   # hollow circle
    "Oaxaca + Chiapas + Upper CA" = 17,  # triangle
    "Lower CA" = 15,       # square
    "Northern Andes" = 18, # diamond
    "Central Andes" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Within-region speciation rate"
  ) +
  theme_classic() + 
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_disp_viburnum_1 = ggplot() +
  geom_line(data = curve_viburnum_df_1,
            aes(x = num_spec, y = total_in_disp_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_1,
            aes(x = num_spec, y = total_in_disp_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_1,
            aes(x = num_spec, y = total_in_disp),
            linewidth = 1) +
  geom_point(data = viburnum_dispersal_dat_1,
             aes(x = num_spec, y = total_in_disp, color = region, shape = region),
             size = 10,
             stroke = 2,
             alpha = 0.8) +
  scale_color_manual(values = c(
    "E/C/W Mexico" = "#d62728",  # red
    "Oaxaca + Chiapas + Upper CA" = "#FFD700",      # yellow
    "Lower CA" = "#17becf",   # cyan
    "Northern Andes" = "#8c2bbd",       # purple
    "Central Andes" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "E/C/W Mexico" = 1,   # hollow circle
    "Oaxaca + Chiapas + Upper CA" = 17,  # triangle
    "Lower CA" = 15,       # square
    "Northern Andes" = 18, # diamond
    "Central Andes" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Total inbound dispersal rate"
  ) +
  theme_classic() + 
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

p_disp_viburnum_2 = ggplot() +
  geom_line(data = curve_viburnum_df_2,
            aes(x = num_spec, y = total_in_disp_low),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_2,
            aes(x = num_spec, y = total_in_disp_high),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(data = curve_viburnum_df_2,
            aes(x = num_spec, y = total_in_disp),
            linewidth = 1) +
  geom_point(data = viburnum_dispersal_dat_2,
             aes(x = num_spec, y = total_in_disp, color = region, shape = region),
             size = 10,
             stroke = 2, 
             alpha = 0.8) +
  scale_color_manual(values = c(
    "E/C/W Mexico" = "#d62728",  # red
    "Oaxaca + Chiapas + Upper CA" = "#FFD700",      # yellow
    "Lower CA" = "#17becf",   # cyan
    "Northern Andes" = "#8c2bbd",       # purple
    "Central Andes" = "#228B22"           # green
  )) + 
  scale_shape_manual(values = c(
    "E/C/W Mexico" = 1,   # hollow circle
    "Oaxaca + Chiapas + Upper CA" = 17,  # triangle
    "Lower CA" = 15,       # square
    "Northern Andes" = 18, # diamond
    "Central Andes" = 3    # plus
  )) +
  labs(
    x = "Num. extant species",
    y = "Total inbound dispersal rate"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    #
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    #
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    #
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 20)
  )+
  scale_y_log10()+
  scale_x_log10()

# save 

pdf(paste0(out_fp,"anolis_emp_plots_network_1.pdf"), width = 10, height = 10)

print(p_spec_anolis_1)
print(p_ext_anolis_1)
print(p_disp_anolis_1)

dev.off()

pdf(paste0(out_fp,"anolis_emp_plots_network_2.pdf"), width = 10, height = 10)

print(p_spec_anolis_2)
print(p_ext_anolis_2)
print(p_disp_anolis_2)

dev.off()

pdf(paste0(out_fp,"viburnum_emp_plots_network_1.pdf"), width = 10, height = 10)

print(p_spec_viburnum_1)
print(p_disp_viburnum_1)

dev.off()

pdf(paste0(out_fp,"viburnum_emp_plots_network_2.pdf"), width = 10, height = 10)

print(p_spec_viburnum_2)
print(p_disp_viburnum_2)

dev.off()
