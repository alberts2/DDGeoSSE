# # This R script plots the average region occupancy across data (trees) in continuous time for 3-region system

library(ggplot2)

plot_reg_occupancy = function(num_sim){
  ###################
  # # SET FILE SYSTEM
  ####################
  fp     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
  region_dat_fp  = paste0(fp,"data/sim1/cts_species_count_reg/") #input region occupancy data in cts time 
  # range_dat_fp  = paste0(fp,"data/sim1/cts_species_count_range/") #input range occupancy data in cts time
  out_fp = paste0(fp,"plot/sim1/") #output 
  
  #
  num_region = 3
  #
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
  # a = function(test){
  for (i in 1:num_sim){
    dat = read.table(paste0(region_dat_fp,"cts_species_count_region_tree_",i,".csv"),header = T)
    
    # create dataframe for that tree
    dat[1,] = gsub("[a-zA-Z]", "", dat[1,])
    columns = lapply(1:num_region,function(i) NA)
    names(columns) <- paste0("Region ", 1:num_region)
    #
    dat_df = data.frame("time"= NA)
    dat_df = cbind(dat_df,as.data.frame(columns))
    
    for (j in 1:nrow(dat)){
      dat[j,] = gsub("[a-zA-Z]", "", dat[j,]) #remove letters
      #
      dat_df[j,] = unlist(strsplit(dat[j,], ",")) #unlist 
      # browser()
      for (k in 1:length(dat_df[j,])){ #turn each element of dat_df[j,] into a numeric value and replace
        if (k == 1){
          dat_df[j,k] = gsub("\\[.*?", "", dat_df[j,k])
        } else if (k == length(dat_df[j,])){
          dat_df[j,k] = gsub("*?\\]", "", dat_df[j,k])
        } else {
          dat_df[j,k] = dat_df[j,k]
        }
      }
    }
    dat_df[, 1:(num_region+1)] <- sapply(dat_df[, 1:(num_region+1)], as.numeric) #convert to numeric
    #
    region_1_arr[[i]] = dat_df$Region.1
    region_2_arr[[i]] = dat_df$Region.2
    region_3_arr[[i]] = dat_df$Region.3
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
  
  #### PLOT 
  dat_reg_1 = data.frame(time = time_vec,mean_reg = mean_reg_1)
  dat_reg_2 = data.frame(time = time_vec,mean_reg = mean_reg_2)
  dat_reg_3 = data.frame(time = time_vec,mean_reg = mean_reg_3)
  #
  p_1 = ggplot(dat_reg_1, aes(x = time_vec, y = mean_reg)) +
        geom_line(color = "blue", size = 1) +
        labs(title = paste0("Average region 1 occupancy across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_2 = ggplot(dat_reg_2, aes(x = time_vec, y = mean_reg)) +
        geom_line(color = "red", size = 1) + 
        labs(title = paste0("Average region 2 occupancy across ",num_sim," trees"), x = "Time", y = "Average species count") 
  
  p_3 = ggplot(dat_reg_3, aes(x = time_vec, y = mean_reg)) +
        geom_line(color = "black", size = 1) +
        labs(title = paste0("Average region 3 occupancy across ",num_sim," trees"), x = "Time", y = "Average species count") 
  plot(p_1)
  plot(p_2)
  plot(p_3)
  # }
  ggsave("region_1_occupancy.jpg",plot = p_1, path = out_fp)
  ggsave("region_2_occupancy.jpg",plot = p_2, path = out_fp)
  ggsave("region_3_occupancy.jpg",plot = p_3, path = out_fp)
  
  return(summary_df)
}
