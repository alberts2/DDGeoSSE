# # This R script plots the average region occupancy across data (trees) in continuous time for 3-region system
# # where the equilibrium diversity across regions can differ 

library(ggplot2)


  ###################
  # # SET FILE SYSTEM
  ####################
  fp     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
  region_dat_fp  = paste0(fp,"data/simulation_areasize/simulate/") #input region occupancy data in cts time 
  out_fp = paste0(fp,"plot/unequal_carrying") #output
  
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
  num_sim = 100
  for (i in 1:num_sim){
    dat_df = read.table(paste0(region_dat_fp,"cts_species_count_region_tree_",i,".csv"),sep = ",",header = T)
    # In the case of tree extinction, fill the rest of entries with 0
    idx_1 = which(is.na(dat_df$region_1))
    idx_2 = which(is.na(dat_df$region_2))
    idx_3 = which(is.na(dat_df$region_3))
    dat_df$region_1[idx_1]=0
    dat_df$region_2[idx_2]=0
    dat_df$region_3[idx_3]=0
    #
    region_1_arr[[i]] = dat_df$region_1
    region_2_arr[[i]] = dat_df$region_2
    region_3_arr[[i]] = dat_df$region_3
    print(paste0("done sim ",i))
  }
  #
  time_vec = dat_df$time
  #
  ci_summary <- function(vals) {
    m  <- mean(vals)
    se <- sd(vals) / sqrt(length(vals))
    
    c(
      mean  = m,
      lower = m - 1.96 * se,
      upper = m + 1.96 * se
    )
  }
  # #
  # for (i in 1:length(time_vec)){
  #   mean_reg_1 = append(mean_reg_1, mean(sapply(region_1_arr,function(x) x[i])))
  #   mean_reg_2 = append(mean_reg_2, mean(sapply(region_2_arr,function(x) x[i])))
  #   mean_reg_3 = append(mean_reg_3, mean(sapply(region_3_arr,function(x) x[i])))
  # }
  # #
  # summary_df = data.frame(time = time_vec,mean_reg_1 = mean_reg_1,
  #                         mean_reg_2 = mean_reg_2, mean_reg_3 = mean_reg_3)
  # 
  # ###################
  # # # PLOT
  # ###################
  # dat_reg_1 = data.frame(time = time_vec,mean_reg = mean_reg_1)
  # dat_reg_2 = data.frame(time = time_vec,mean_reg = mean_reg_2)
  # dat_reg_3 = data.frame(time = time_vec,mean_reg = mean_reg_3)
  sum_reg_1 <- t(sapply(seq_along(time_vec), function(i) {
    vals <- sapply(region_1_arr, function(x) x[i])
    ci_summary(vals)
  }))
  
  sum_reg_2 <- t(sapply(seq_along(time_vec), function(i) {
    vals <- sapply(region_2_arr, function(x) x[i])
    ci_summary(vals)
  }))
  
  sum_reg_3 <- t(sapply(seq_along(time_vec), function(i) {
    vals <- sapply(region_3_arr, function(x) x[i])
    ci_summary(vals)
  }))
  
  dat_reg_1 <- data.frame(time = time_vec, sum_reg_1)
  dat_reg_2 <- data.frame(time = time_vec, sum_reg_2)
  dat_reg_3 <- data.frame(time = time_vec, sum_reg_3)
  #
  p_1 <- ggplot(dat_reg_1, aes(x = time, y = mean)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2, fill = "blue") +
    geom_line(color = "blue", size = 1) +
    # geom_hline(yintercept = 5,
    #            linetype = "dashed", color = "purple", size = 1) +
    labs(title = paste0("Average region 1 occupancy across ", num_sim, " trees"),
         x = "Time", y = "Average species count") +
    coord_cartesian(xlim = c(0, 500000),ylim = c(0, 25)) +
    theme_minimal()+
    theme(
      axis.text  = element_text(size = 23),
      axis.title = element_text(size = 26),
      plot.title = element_text(size = 26)
    )
  
  p_2 <- ggplot(dat_reg_2, aes(x = time, y = mean)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2, fill = "red") +
    geom_line(color = "red", size = 1) +
    # geom_hline(yintercept = 10,
    #            linetype = "dashed", color = "purple", size = 1) +
    labs(title = paste0("Average region 2 occupancy across ", num_sim, " trees"),
         x = "Time", y = "Average species count") +
    coord_cartesian(xlim = c(0, 500000),ylim = c(0, 25)) +
    theme_minimal()+
    theme(
      axis.text  = element_text(size = 23),
      axis.title = element_text(size = 26),
      plot.title = element_text(size = 26)
    )
  
  p_3 <- ggplot(dat_reg_3, aes(x = time, y = mean)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2, fill = "black") +
    geom_line(color = "black", size = 1) +
    # geom_hline(yintercept = 15,
    #            linetype = "dashed", color = "purple", size = 1) +
    labs(title = paste0("Average region 3 occupancy across ", num_sim, " trees"),
         x = "Time", y = "Average species count") +
    coord_cartesian(xlim = c(0, 500000),ylim = c(0, 25)) +
    theme_minimal()+
    theme(
      axis.text  = element_text(size = 23),
      axis.title = element_text(size = 26),
      plot.title = element_text(size = 26)
    )
  
  
  plot(p_1)
  plot(p_2)
  plot(p_3)
  
  ###################
  # # SAVE PLOT
  ###################
  
  pdf(file = paste0(out_fp,"/region_occupancy_unequal.pdf"), width = 9, height = 6)
  print(p_1)
  print(p_2)
  print(p_3)
  dev.off()
  
