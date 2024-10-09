# # This R script plots the average region occupancy across data (trees) in continuous time for 3-region system

# # TO DO: No need to post-process the csv anymore, need to edit this script

library(ggplot2)

plot_range_occupancy = function(num_sim,batch_num){
  ###################
  # # SET FILE SYSTEM
  ###################
  fp     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
  range_dat_fp  = paste0(fp,"data/sim2/cts_species_count_range/") #input range occupancy data in cts time
  out_fp = paste0(fp,"plot/sim2/") #output 
  
  ###################
  # # INITIALIZATION
  ###################
  range_A_arr = list() # list containing range A occupancy in cts time for each tree 
  range_B_arr = list() # list containing range B occupancy in cts time for each tree
  range_C_arr = list() # list containing range C occupancy in cts time for each tree
  range_AB_arr = list() # list containing range AB occupancy in cts time for each tree
  range_AC_arr = list() # list containing range AC occupancy in cts time for each tree
  range_BC_arr = list() # list containing range BC occupancy in cts time for each tree
  range_ABC_arr = list() # list containing range ABC occupancy in cts time for each tree
  
  mean_range_1   = c() # vector containing average range A occupancy at each time point across trees 
  mean_range_2   = c() # vector containing average range B occupancy at each time point across trees 
  mean_range_3   = c() # vector containing average range C occupancy at each time point across trees
  mean_range_4   = c() # vector containing average range AB occupancy at each time point across trees 
  mean_range_5   = c() # vector containing average range AC occupancy at each time point across trees 
  mean_range_6   = c() # vector containing average range BC occupancy at each time point across trees 
  mean_range_7   = c() # vector containing average range ABC occupancy at each time point across trees 
  
  ###################
  # # READ DATA
  ###################
  num_sim = num_sim
  # a = function(test){
  for (i in 1:num_sim){
    dat_df = read.table(paste0(range_dat_fp,"cts_species_count_range_tree_",i,".csv"),sep = ",",header = T)
    #
    range_A_arr[[i]] = dat_df$range_1
    range_B_arr[[i]] = dat_df$range_2
    range_C_arr[[i]] = dat_df$range_3
    range_AB_arr[[i]] = dat_df$range_4
    range_AC_arr[[i]] = dat_df$range_5
    range_BC_arr[[i]] = dat_df$range_6
    range_ABC_arr[[i]] = dat_df$range_7
    print(paste0("done sim ",i))
  }
  #
  time_vec = dat_df$time
  #
  for (i in 1:length(time_vec)){
    mean_range_1 = append(mean_range_1, mean(sapply(range_A_arr,function(x) x[i])))
    mean_range_2 = append(mean_range_2, mean(sapply(range_B_arr,function(x) x[i])))
    mean_range_3 = append(mean_range_3, mean(sapply(range_C_arr,function(x) x[i])))
    mean_range_4 = append(mean_range_4, mean(sapply(range_AB_arr,function(x) x[i])))
    mean_range_5 = append(mean_range_5, mean(sapply(range_AC_arr,function(x) x[i])))
    mean_range_6 = append(mean_range_6, mean(sapply(range_BC_arr,function(x) x[i])))
    mean_range_7 = append(mean_range_7, mean(sapply(range_ABC_arr,function(x) x[i])))
    
  }
  #
  summary_df = data.frame(time = time_vec,mean_range_1 = mean_range_1,mean_range_2 = mean_range_2,
                          mean_range_3 = mean_range_3,mean_range_4 = mean_range_4,mean_range_5 = mean_range_5,
                          mean_range_6 = mean_range_6,mean_range_7 = mean_range_7)
  
  ###################
  # # PLOT
  ###################
  dat_range_1 = data.frame(time = time_vec,mean_range = mean_range_1)
  dat_range_2 = data.frame(time = time_vec,mean_range = mean_range_2)
  dat_range_3 = data.frame(time = time_vec,mean_range = mean_range_3)
  dat_range_4 = data.frame(time = time_vec,mean_range = mean_range_4)
  dat_range_5 = data.frame(time = time_vec,mean_range = mean_range_5)
  dat_range_6 = data.frame(time = time_vec,mean_range = mean_range_6)
  dat_range_7 = data.frame(time = time_vec,mean_range = mean_range_7)
  #
  p_1 = ggplot(dat_range_1, aes(x = time_vec, y = mean_range)) +
    geom_line(color = "blue", size = 1) +
    labs(title = paste0("Average species count with range A across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_2 = ggplot(dat_range_2, aes(x = time_vec, y = mean_range)) +
    geom_line(color = "red", size = 1) +
    labs(title = paste0("Average species count with range B across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_3 = ggplot(dat_range_3, aes(x = time_vec, y = mean_range)) +
    geom_line(color = "black", size = 1) +
    labs(title = paste0("Average species count with range C across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_4 = ggplot(dat_range_4, aes(x = time_vec, y = mean_range)) +
    geom_line(color = "orange", size = 1) +
    labs(title = paste0("Average species count with range AB across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_5 = ggplot(dat_range_5, aes(x = time_vec, y = mean_range)) +
    geom_line(color = "green", size = 1) +
    labs(title = paste0("Average species count with range AC across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_6 = ggplot(dat_range_6, aes(x = time_vec, y = mean_range)) +
    geom_line(color = "yellow", size = 1) +
    labs(title = paste0("Average species count with range BC across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_7 = ggplot(dat_range_7, aes(x = time_vec, y = mean_range)) +
    geom_line(color = "brown", size = 1) +
    labs(title = paste0("Average species count with range ABC across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  plot(p_1)
  plot(p_2)
  plot(p_3)
  plot(p_4)
  plot(p_5)
  plot(p_6)
  plot(p_7)
  
  ###################
  # # SAVE PLOTS
  ###################
  ggsave("range_A_count.jpg",plot = p_1, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("range_B_count.jpg",plot = p_2, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("range_C_count.jpg",plot = p_3, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("range_AB_count.jpg",plot = p_4, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("range_AC_count.jpg",plot = p_5, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("range_BC_count.jpg",plot = p_6, path = out_fp,width = 8, height = 6, units = "in")
  ggsave("range_ABC_count.jpg",plot = p_7, path = out_fp,width = 8, height = 6, units = "in")
  # }
  return(summary_df)
}
