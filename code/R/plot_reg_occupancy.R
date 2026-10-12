# # This R script plots the average region occupancy across data (trees) in continuous time for 3-region system

library(ggplot2)

plot_reg_occupancy = function(num_sim,lemma){
  ###################
  # # SET FILE SYSTEM
  ####################
  fp     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
  region_dat_fp  = paste0(fp,"data/theory_validation/lemma",lemma,"_new/cts_species_count_reg/") #input region occupancy data in cts time 
  out_fp = paste0(fp,"plot/lemma",lemma,"_new/") #output
  
  ###################
  # # INITIALIZATION
  ####################
  region_1_arr = list() # list containing region 1 occupancy in cts time for each tree 
  region_2_arr = list() # list containing region 2 occupancy in cts time for each tree 
  region_3_arr = list() # list containing region 1 occupancy in cts time for each tree 
  mean_reg_1   = c() # vector containing average region 1 occupancy at each time point across trees 
  mean_reg_2   = c() # vector containing average region 2 occupancy at each time point across trees 
  mean_reg_3   = c() # vector containing average region 3 occupancy at each time point across trees 
  
  ###################
  # # READ DATA
  ####################
  num_sim = num_sim
  for (i in 1:num_sim){
    dat_df = read.table(paste0(region_dat_fp,"cts_species_count_region_tree_",i,".csv"),sep = ",",header = T)
    #
    region_1_arr[[i]] = dat_df$region_1
    region_2_arr[[i]] = dat_df$region_2
    region_3_arr[[i]] = dat_df$region_3
    print(paste0("done sim ",i))
  }
  #
  time_vec = dat_df$time
  #
  for (i in 1:length(time_vec)){
    mean_reg_1 = append(mean_reg_1, mean(sapply(region_1_arr,function(x) x[i])))
    mean_reg_2 = append(mean_reg_2, mean(sapply(region_2_arr,function(x) x[i])))
    mean_reg_3 = append(mean_reg_3, mean(sapply(region_3_arr,function(x) x[i])))
  }
  #
  summary_df = data.frame(time = time_vec,mean_reg_1 = mean_reg_1,
                          mean_reg_2 = mean_reg_2, mean_reg_3 = mean_reg_3)
  
  ###################
  # # PLOT
  ###################
  dat_reg_1 = data.frame(time = time_vec,mean_reg = mean_reg_1)
  dat_reg_2 = data.frame(time = time_vec,mean_reg = mean_reg_2)
  dat_reg_3 = data.frame(time = time_vec,mean_reg = mean_reg_3)
  #
  p_1 = ggplot(dat_reg_1, aes(x = time_vec, y = mean_reg)) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = 19, linetype = "dashed", color = "purple", size = 1)+
    labs(title = paste0("Average region 1 occupancy across ",num_sim," trees"), x = "Time", y = "Average species count")+
    theme_minimal() +theme(axis.text = element_text(size = 26))
  
  p_2 = ggplot(dat_reg_2, aes(x = time_vec, y = mean_reg)) +
    geom_line(color = "red", size = 1) + 
    geom_hline(yintercept = 19, linetype = "dashed", color = "purple", size = 1)+
    labs(title = paste0("Average region 2 occupancy across ",num_sim," trees"), x = "Time", y = "Average species count")+
    theme_minimal() +theme(axis.text = element_text(size = 26))
  
  p_3 = ggplot(dat_reg_3, aes(x = time_vec, y = mean_reg)) +
    geom_line(color = "black", size = 1) +
    geom_hline(yintercept = 19, linetype = "dashed", color = "purple", size = 1)+
    labs(title = paste0("Average region 3 occupancy across ",num_sim," trees"), x = "Time", y = "Average species count")+
    theme_minimal() +theme(axis.text = element_text(size = 26))
  
  plot(p_1)
  plot(p_2)
  plot(p_3)
  
  ###################
  # # SAVE PLOT
  ###################
  ggsave("region_1_occupancy.jpg",plot = p_1, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("region_2_occupancy.jpg",plot = p_2, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("region_3_occupancy.jpg",plot = p_3, path = out_fp,width = 8, height = 6, units = "in")
  # 
  return(summary_df)
}
