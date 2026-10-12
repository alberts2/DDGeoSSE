# # This R script plots the average region occupancy across data (trees) in continuous time for 3-region system

library(ggplot2)
library(gridExtra)
# plot_conceptual = function(num_sim){
###################
# # SET FILE SYSTEM
####################
fp     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
region_dat_fp  = paste0(fp,"data/conceptual_fig/") #input region occupancy data in cts time
# region_dat_fp  = paste0(fp,"data/test/") #input region occupancy data in cts time
out_fp = paste0(fp,"plot/conceptual_fig/") #output 

###################
# # INITIALIZATION
####################
region_1_arr = list() # list containing region 1 occupancy in cts time for each tree 
region_2_arr = list() # list containing region 2 occupancy in cts time for each tree 
mean_reg_1   = c() # vector containing average region 1 occupancy at each time point across trees 
mean_reg_2   = c() # vector containing average region 2 occupancy at each time point across trees 


###################
# # READ DATA
####################
num_sim = 1
num_sim = num_sim
for (i in 1:num_sim){
  j = i - 1
  dat_df = read.table(paste0(region_dat_fp,"species_count_region_tree_",j,".csv"),sep = ",",header = T)
  # dat_df = read.table(paste0(region_dat_fp,"cts_species_count_region_tree_",j,".csv"),sep = ",",header = T)
  #
  region_1_arr[[i]] = dat_df$region_1
  region_2_arr[[i]] = dat_df$region_2
  print(paste0("done sim ",i))
}
#
time_vec = dat_df$time
#
for (i in 1:length(time_vec)){
  mean_reg_1 = append(mean_reg_1, mean(sapply(region_1_arr,function(x) x[i])))
  mean_reg_2 = append(mean_reg_2, mean(sapply(region_2_arr,function(x) x[i])))
}
#
# total number of extant species through time 
total_species = mean_reg_1 + mean_reg_2 
# read the simulating parameter 
param_df = read.table(paste0(region_dat_fp,"sim.",j,".labels.csv"),sep = ",",header = T)
# list the params
rho_w = exp(param_df$log_base_w)
rho_e = exp(param_df$log_base_e)
rho_d = exp(param_df$log_base_d)
w_d = log(param_df$exp_scalar_w)
e_d = log(param_df$exp_scalar_e)
d_d_home = log(param_df$exp_scalar_d_i)
d_d_away = log(param_df$exp_scalar_d_j)
# 
spec_rate_1_vec = rep(NA,length(mean_reg_1))
spec_rate_2_vec = rep(NA,length(mean_reg_1))
ext_rate_1_vec = rep(NA,length(mean_reg_1))
ext_rate_2_vec = rep(NA,length(mean_reg_1))
disp_rate_1_2_vec = rep(NA,length(mean_reg_1))
disp_rate_2_1_vec = rep(NA,length(mean_reg_1))
#
for (i in 1:length(mean_reg_1)){
  if (mean_reg_1[i] == 0){
    spec_rate_1_vec[i] = 0
    ext_rate_1_vec[i] = 0
    print("1")
    #
    spec_rate_2_vec[i] = rho_w*(log(mean_reg_2[i])+1)^{w_d}
    ext_rate_2_vec[i] = rho_e*(log(mean_reg_2[i])+1)^{e_d}
  } else if (mean_reg_2[i] == 0){
    spec_rate_1_vec[i] = rho_w*(log(mean_reg_1[i])+1)^{w_d}
    ext_rate_1_vec[i] = rho_e*(log(mean_reg_1[i])+1)^{e_d}
    #
    spec_rate_2_vec[i] = 0
    ext_rate_2_vec[i] = 0
    print("2")
  } else if (mean_reg_1[i] != 0 && mean_reg_2[i] != 0){
    spec_rate_1_vec[i] = rho_w*(log(mean_reg_1[i])+1)^{w_d}
    spec_rate_2_vec[i] = rho_w*(log(mean_reg_2[i])+1)^{w_d}
    ext_rate_1_vec[i] = rho_e*(log(mean_reg_1[i])+1)^{e_d}
    ext_rate_2_vec[i] = rho_e*(log(mean_reg_2[i])+1)^{e_d}
  }
  if (mean_reg_1[i]==0){
    disp_rate_1_2_vec[i] = 0
  } else {
    disp_rate_1_2_vec[i] = rho_d*(log(mean_reg_1[i])+1)^{d_d_home}*(log(mean_reg_2[i]+1)+1)^{d_d_away}
  }
  if (mean_reg_2[i]==0){
    disp_rate_2_1_vec[i] = 0
  } else {
    disp_rate_2_1_vec[i] = rho_d*(log(mean_reg_2[i])+1)^{d_d_home}*(log(mean_reg_1[i]+1)+1)^{d_d_away}
  }
}
#
summary_df = data.frame(time = time_vec,mean_reg_1 = mean_reg_1,
                        mean_reg_2 = mean_reg_2, tot_spec = total_species,
                        spec_rate_1 = spec_rate_1_vec, spec_rate_2 = spec_rate_2_vec,
                        ext_rate_1 = ext_rate_1_vec, ext_rate_2 = ext_rate_2_vec,
                        disp_rate_1_2 = disp_rate_1_2_vec, disp_rate_2_1 = disp_rate_2_1_vec)

# Replace all 0 rates with NA 
summary_df$spec_rate_1[summary_df$spec_rate_1==0]     = NA
summary_df$spec_rate_2[summary_df$spec_rate_2==0]     = NA
summary_df$ext_rate_1[summary_df$ext_rate_1==0]       = NA
summary_df$ext_rate_2[summary_df$ext_rate_2==0]       = NA
summary_df$disp_rate_1_2[summary_df$disp_rate_1_2==0] = NA
summary_df$disp_rate_2_1[summary_df$disp_rate_2_1==0] = NA

# add symbol at the first non-zero rate from each process
val_spec_1 = summary_df[which(!is.na(summary_df$spec_rate_1))[1], ]
val_spec_2 = summary_df[which(!is.na(summary_df$spec_rate_2))[1], ]
val_ext_1  = summary_df[which(!is.na(summary_df$ext_rate_1))[1], ]
val_ext_2  = summary_df[which(!is.na(summary_df$ext_rate_2))[1], ]
val_disp_1 = summary_df[which(!is.na(summary_df$disp_rate_1_2))[1], ]
val_disp_2 = summary_df[which(!is.na(summary_df$disp_rate_2_1))[1], ]
###################
# # PLOT
###################
dat_reg_1 = data.frame(time = time_vec,num_spec = mean_reg_1)
dat_reg_2 = data.frame(time = time_vec,num_spec = mean_reg_2)
# Plot region occupancy
p_1 = ggplot(summary_df, aes(x = time_vec)) +
  geom_line(aes(y = mean_reg_1, color = "species_1"), size = 1) +
  geom_line(aes(y = mean_reg_2, color = "species_2"), size = 1) + 
  geom_line(aes(y = tot_spec, color = "tot_species"), size = 1,alpha = 0.8) + 
  scale_y_continuous(limits = c(0, 10),breaks = seq(0, 10, by = 2))+
  scale_color_manual(values = c("species_1" = "blue", "species_2" = "red", "tot_species" = "gray")) +
  labs(title = paste0("Region occupancy"), x = "Time", y = "Species count") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 20),       # tick labels
    legend.position = "none"
  )

# Plot speciation rates in both regions
p_2 = ggplot(summary_df, aes(x = time)) +
  geom_line(aes(y = spec_rate_1, color = "Region 1"), size = 1) +
  geom_line(aes(y = spec_rate_2, color = "Region 2"), size = 1) +
  geom_point(
    data = val_spec_1,
    aes(x = time,y = spec_rate_1),
    shape = 8,
    color = "blue",
    size = 4
  ) +
  geom_point(
    data = val_spec_2,
    aes(x = time,y = spec_rate_2),
    shape = 8,
    color = "red",
    size = 4
  ) +
  scale_y_log10(limits = c(0.01, 1)) +
  scale_color_manual(values = c("Region 1" = "blue", "Region 2" = "red")) +
  labs(title = "Speciation rate in both regions",
       x = "Time",
       y = "Speciation rate",
       color = "Region") + 
  theme_minimal()+
  theme(
    axis.text = element_text(size = 20),       # tick labels
    legend.position = "none"
  ) 

# Plot extinction rates in both regions
p_3 = ggplot(summary_df, aes(x = time)) +
  geom_line(aes( y = ext_rate_1, color = "Region 1"), size = 1) + 
  geom_line(aes( y = ext_rate_2, color = "Region 2"), size = 1) + 
  geom_point(
    data = val_ext_1,
    aes(x = time,y = ext_rate_1),
    shape = 8,
    color = "blue",
    size = 4
  ) +
  geom_point(
    data = val_ext_2,
    aes(x = time,y = ext_rate_2),
    shape = 8,
    color = "red",
    size = 4
  ) +
  scale_y_log10(limits = c(0.01, 1))+
  scale_color_manual(values = c("Region 1" = "blue", "Region 2" = "red")) +
  labs(title = "Extinction rate in both regions",
       x = "Time",
       y = "Extinction rate",
       color = "Region") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 20),       # tick labels
    legend.position = "none"
  ) 

# Plot dispersal rates in both regions
p_4 = ggplot(summary_df, aes(x = time)) +
  geom_line(aes(y = disp_rate_2_1, color = "Into 1"), size = 1) + 
  geom_line(aes(y = disp_rate_1_2, color = "Into 2"), size = 1) + 
  geom_point(
    data = val_disp_2,
    aes(x = time,y = disp_rate_2_1),
    shape = 8,
    color = "blue",
    size = 4
  ) +
  geom_point(
    data = val_disp_1,
    aes(x = time,y = disp_rate_1_2),
    shape = 8,
    color = "red",
    size = 4
  ) +
  scale_y_log10(limits = c(0.01, 1))+
  scale_color_manual(values = c("Into 1" = "blue", "Into 2" = "red")) +
  labs(title = "Dispersal rate into 1 and into 2",
       x = "Time",
       y = "Dispersal rate",
       color = "Region") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 20),       # tick labels
    legend.position = "none"
  ) 
 

p_5 = grid.arrange(p_1, p_2, p_3, p_4, ncol = 2)

###################
# # SAVE PLOT
###################
ggsave("region_occupancy.jpg",plot = p_1, path = out_fp,width = 8, height = 4, units = "in")
ggsave("speciation_rate.jpg",plot = p_2, path = out_fp,width = 8, height = 4, units = "in")
ggsave("extinction_rate.jpg",plot = p_3, path = out_fp,width = 8, height = 4, units = "in")
ggsave("dispersal_rate.jpg",plot = p_4, path = out_fp,width = 8, height = 4, units = "in")
ggsave("combined_plot.jpg",plot = p_5, path = out_fp,width = 8, height = 4, units = "in")
# return(summary_df)
# }

