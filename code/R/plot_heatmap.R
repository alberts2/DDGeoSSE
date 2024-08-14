library(ggplot2)

##############
# EXTINCTION  #
##############

# LEGEND 
# rate_scale = rate scaling for the desity component on extinction process (-inf,inf)
# base rate  = base rate on (rho) for extinction (rho > 0)
# indep_rate = density-independent rate for extinction (> 0)
# n1_lower   = mininum number of species in region i (> 0)
# n2_lower   = minimum number of species in region j (> 0)
# timewin    = number of timebins

# test = heat_map_e(rate_scale = 0.1, base_rate = 0.1, indep_rate_vec = c(0.1,0.1),
#                  n1_lower = 10, n2_lower = 20, timewin = 20)

heat_map_e <- function(rate_scale,base_rate,indep_rate_vec,n1_lower,n2_lower,timewin){
  ##############
  # CONDITION CHECK  #
  ##############
  #
  if(n1_lower < 0 || n2_lower < 0){
    error("Number of species in any region cannot be negative")
  }
  # browser()
  #
  n1_vec     <- seq(n1_lower,n1_lower+timewin-1,by=1)
  n2_vec     <- seq(n2_lower,n2_lower+timewin-1,by=1)
  # create grid 
  n_grid <- expand.grid(n1 = n1_vec, n2 = n2_vec)
  # browser()
  #
  # append density dependent effect and rate over time 
  den_vec_log     <- vector("list",2) # log density effect
  den_vec_loglog  <- vector("list",2) # loglog density effect
  #
  ei_vec_log     <- vector("list",2)  # log density-dependent extinction
  ei_vec_loglog  <- vector("list",2)  # loglog density-dependent extinction
  #
  # browser()
  for (i in 1:2){ # loop over the two regions
    for (j in 1:length(n_grid$n1)){ # loop over rows in n_grid
      den_vec_log[[i]][j]     <- n_grid[j,i]^rate_scale 
      den_vec_loglog[[i]][j]  <- exp(rate_scale*log(log(n_grid[j,i])+1))
      #
      ei_vec_log[[i]][j]     <- base_rate*indep_rate_vec[i] * den_vec_log[[i]][j]
      ei_vec_loglog[[i]][j]  <- base_rate*indep_rate_vec[i] * den_vec_loglog[[i]][j]
    }
  }
  # browser()
  ##############
  # DATA PROCESSING   #
  ##############
  dat_log             <- data.frame(func = rep("log",length(n_grid$n1)),num_reg_1 = n_grid[,1], num_reg_2 = n_grid[,2],
                                    rate_scale = rep(rate_scale,length(n_grid$n1)),
                                    indep_rate_1 = rep(indep_rate_vec[1],length(n_grid$n1)),  indep_rate_2 = rep(indep_rate_vec[2],length(n_grid$n1)),
                                    density_effect_1 = den_vec_log[[1]],
                                    density_effect_2 = den_vec_log[[2]],
                                    overall_rate_1 = ei_vec_log[[1]],overall_rate_2 = ei_vec_log[[2]],
                                    rate_ratio_1_2 = ei_vec_log[[1]]/ei_vec_log[[2]])
  

  # browser()
  dat_loglog          <- data.frame(func = rep("loglog",length(n_grid$n1)),num_reg_1 = n_grid[,1], num_reg_2 = n_grid[,2],
                                    rate_scale = rep(rate_scale,length(n_grid$n1)),
                                    indep_rate_1 = rep(indep_rate_vec[1],length(n_grid$n1)),  indep_rate_2 = rep(indep_rate_vec[2],length(n_grid$n1)),
                                    density_effect_1 = den_vec_loglog[[1]],
                                    density_effect_2 = den_vec_loglog[[2]],
                                    overall_rate_1 = ei_vec_loglog[[1]],overall_rate_2 = ei_vec_loglog[[2]],
                                    rate_ratio_1_2 = ei_vec_loglog[[1]]/ei_vec_loglog[[2]])
  
  #
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog)
  #
  #
  ##############
  # PLOTTINGS  #
  ##############
  func_names <- c("log","loglog")
  plot_list = list()
  title_list <- c(paste("log extinction\n","rate scale=",rate_scale,",","base rate=",base_rate,"\n",
                        "independent rate 1=",indep_rate_vec[1],",","independent rate 2=",indep_rate_vec[2]),
                  paste("log-log extinction\n","rate scale=",rate_scale,",","base rate=",base_rate,"\n",
                        "independent rate 1=",indep_rate_vec[1],",","independent rate 2=",indep_rate_vec[2]))
  # make plot for each model
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    # browser()
    # browser()
    #
    p = ggplot(dat_i, aes(x=num_reg_1 ,y=num_reg_2,fill=rate_ratio_1_2))
    #
    p = p+geom_tile()
    # plot labels 
    p = p + xlab("#species in region 1") + ylab("#species in region 2") + ggtitle(title_list[i]) + theme(plot.title = element_text(hjust = 0.5))
    #
    plot_list[[i]] = p
    print(p)
  }
  plot_list[[1]]
  plot_list[[2]]
  #
  return(dat_comb)
}

##############
# WITHIN REGION SPECIATION #
##############

# LEGEND 
# rate_scale = rate scaling for the desity component on within-reg speciation process (-inf,inf)
# base rate  = base rate on (rho) for within-reg speciation process (rho > 0)
# indep_rate = density-independent rate for extinction (> 0)
# n1_lower   = mininum number of species in region i (> 0)
# n2_lower   = minimum number of species in region j (> 0)
# timewin    = number of timebins

heat_map_w <- function(rate_scale,base_rate,indep_rate_vec,n1_lower,n2_lower,timewin){
  ##############
  # CONDITION CHECK  #
  ##############
  #
  if(n1_lower < 0 || n2_lower < 0){
    error("Number of species in any region cannot be negative")
  }
  # browser()
  #
  n1_vec     <- seq(n1_lower,n1_lower+timewin-1,by=1)
  n2_vec     <- seq(n2_lower,n2_lower+timewin-1,by=1)
  # create grid 
  n_grid <- expand.grid(n1 = n1_vec, n2 = n2_vec)
  # browser()
  #
  # append density dependent effect and rate over time 
  den_vec_log     <- vector("list",2) # log density effect
  den_vec_loglog  <- vector("list",2) # loglog density effect
  #
  wi_vec_log     <- vector("list",2)  # log density-dependent within-region speciation
  wi_vec_loglog  <- vector("list",2)  # loglog density-dependent within-region speciation
  #
  # browser()
  for (i in 1:2){ # loop over the two regions
    for (j in 1:length(n_grid$n1)){ # loop over rows in n_grid
      den_vec_log[[i]][j]     <- n_grid[j,i]^rate_scale 
      den_vec_loglog[[i]][j]  <- exp(rate_scale*log(log(n_grid[j,i])+1))
      #
      wi_vec_log[[i]][j]     <- base_rate*indep_rate_vec[i] * den_vec_log[[i]][j]
      wi_vec_loglog[[i]][j]  <- base_rate*indep_rate_vec[i] * den_vec_loglog[[i]][j]
    }
  }
  # browser()
  ##############
  # DATA PROCESSING   #
  ##############
  dat_log             <- data.frame(func = rep("log",length(n_grid$n1)),num_reg_1 = n_grid[,1], num_reg_2 = n_grid[,2],
                                    rate_scale = rep(rate_scale,length(n_grid$n1)),
                                    indep_rate_1 = rep(indep_rate_vec[1],length(n_grid$n1)),  indep_rate_2 = rep(indep_rate_vec[2],length(n_grid$n1)),
                                    density_effect_1 = den_vec_log[[1]],
                                    density_effect_2 = den_vec_log[[2]],
                                    overall_rate_1 = wi_vec_log[[1]],overall_rate_2 = wi_vec_log[[2]],
                                    rate_ratio_1_2 = wi_vec_log[[1]]/wi_vec_log[[2]])
  
  
  # browser()
  dat_loglog          <- data.frame(func = rep("loglog",length(n_grid$n1)),num_reg_1 = n_grid[,1], num_reg_2 = n_grid[,2],
                                    rate_scale = rep(rate_scale,length(n_grid$n1)),
                                    indep_rate_1 = rep(indep_rate_vec[1],length(n_grid$n1)),  indep_rate_2 = rep(indep_rate_vec[2],length(n_grid$n1)),
                                    density_effect_1 = den_vec_loglog[[1]],
                                    density_effect_2 = den_vec_loglog[[2]],
                                    overall_rate_1 = wi_vec_loglog[[1]],overall_rate_2 = wi_vec_loglog[[2]],
                                    rate_ratio_1_2 = wi_vec_loglog[[1]]/wi_vec_loglog[[2]])
  
  #
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog)
  #
  #
  ##############
  # PLOTTINGS  #
  ##############
  func_names <- c("log","loglog")
  plot_list = list()
  title_list <- c(paste("log within-region speciation\n","rate scale=",rate_scale,",","base rate=",base_rate,"\n",
                        "independent rate 1=",indep_rate_vec[1],",","independent rate 2=",indep_rate_vec[2]),
                  paste("log-log within-region speciation\n","rate scale=",rate_scale,",","base rate=",base_rate,"\n",
                        "independent rate 1=",indep_rate_vec[1],",","independent rate 2=",indep_rate_vec[2]))
  # make plot for each model
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    # browser()
    # browser()
    #
    p = ggplot(dat_i, aes(x=num_reg_1 ,y=num_reg_2,fill=rate_ratio_1_2))
    #
    p = p+geom_tile()
    # plot labels 
    p = p + xlab("#species in region 1") + ylab("#species in region 2") + ggtitle(title_list[i]) + theme(plot.title = element_text(hjust = 0.5))
    #
    plot_list[[i]] = p
    print(p)
  }
  plot_list[[1]]
  plot_list[[2]]
  #
  return(dat_comb)
}

##############
# DISPERSAL #
##############

# LEGEND 
# rate_scale_i = rate scaling for the desity component on source region i (-inf,inf)
# rate_scale_i = rate scaling for the desity component on target region j (-inf,inf)
# base rate  = base rate on (rho) for dispersal (rho > 0)
# indep_rate = density-independent rate for dispersal (> 0)
# n1_lower   = mininum number of species in region i (> 0)
# n2_lower   = minimum number of species in region j (> 0)
# timewin    = number of timebins

heat_map_d <- function(rate_scale_i,rate_scale_j,base_rate,indep_rate,n1_lower,n2_lower,timewin){
  ##############
  # CONDITION CHECK  #
  ##############
  #
  if(n1_lower < 0 || n2_lower < 0){
    error("Number of species in any region cannot be negative")
  }
  # browser()
  #
  n1_vec     <- seq(n1_lower,n1_lower+timewin-1,by=1)
  n2_vec     <- seq(n2_lower,n2_lower+timewin-1,by=1)
  # create grid 
  n_grid <- expand.grid(n1 = n1_vec, n2 = n2_vec)
  # browser()
  #
  # append density dependent effect and rate over time 
  den_vec_log     <- rep(NA,length(n_grid$n1)) # log density effect
  den_vec_loglog  <- rep(NA,length(n_grid$n1)) # loglog density effect
  #
  dij_vec_log     <- rep(NA,length(n_grid$n1))  # log density-dependent dispersal rate i to j
  dij_vec_loglog  <- rep(NA,length(n_grid$n1))  # loglog density-dependent dispersal rate i to j
  #
  # browser()
    for (i in 1:length(n_grid$n1)){ # loop over rows in n_grid
      den_vec_log[i]          <- exp(rate_scale_i*log(n_grid[i,1])+rate_scale_j*log(n_grid[i,2]+1))
      den_vec_loglog[i]       <- exp(rate_scale_i*log(log(n_grid[i,1])+1)+rate_scale_j*log(log(n_grid[i,2]+1)+1))
      #
      dij_vec_log[i]          <- base_rate*indep_rate * den_vec_log[i]
      dij_vec_loglog[i]       <- base_rate*indep_rate * den_vec_loglog[i]
    }
  # browser()
  ##############
  # DATA PROCESSING   #
  ##############
  dat_log             <- data.frame(func = rep("log",length(n_grid$n1)),num_reg_1 = n_grid[,1], num_reg_2 = n_grid[,2],
                                    rate_scale_i = rep(rate_scale_i,length(n_grid$n1)), rate_scale_j = rep(rate_scale_j,length(n_grid$n1)),
                                    indep_rate = rep(indep_rate,length(n_grid$n1)),
                                    density_effect  = den_vec_log,
                                    overall_rate = dij_vec_log)
  
  
  # browser()
  dat_loglog          <- data.frame(func = rep("loglog",length(n_grid$n1)),num_reg_1 = n_grid[,1], num_reg_2 = n_grid[,2],
                                    rate_scale_i = rep(rate_scale_i,length(n_grid$n1)), rate_scale_j = rep(rate_scale_j,length(n_grid$n1)),
                                    indep_rate = rep(indep_rate,length(n_grid$n1)),
                                    density_effect  = den_vec_loglog,
                                    overall_rate = dij_vec_loglog)
  
  #
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog)
  #
  #
  ##############
  # PLOTTINGS  #
  ##############
  func_names <- c("log","loglog")
  plot_list = list()
  title_list <- c(paste("log dispersal rate\n","rate scale i=",rate_scale_i,",","rate scale j=",rate_scale_j,"\n",
                        "base rate=",base_rate,",","independent rate=",indep_rate),
                  paste("log-log dispersal rate\n","rate scale i=",rate_scale_i,",","rate scale j=",rate_scale_j,"\n",
                        "base rate=",base_rate,",","independent rate=",indep_rate))
  # make plot for each model
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    # browser()
    # browser()
    #
    p = ggplot(dat_i, aes(x=num_reg_1 ,y=num_reg_2,fill=overall_rate))
    #
    p = p+geom_tile()
    # plot labels 
    p = p + xlab("#species in region 1") + ylab("#species in region 2") + ggtitle(title_list[i]) + theme(plot.title = element_text(hjust = 0.5))
    #
    plot_list[[i]] = p
    print(p)
  }
  plot_list[[1]]
  plot_list[[2]]
  #
  return(dat_comb)
}

##############
# BETWEEN-REGION SPECIATION #
##############

# LEGEND 
# rate_scale   = rate scaling for the desity component on between-reg speciation process (-inf,inf)
# base rate    = base rate on (rho) for between-reg speciation process (rho > 0)
# left_range   = left daughter species' range 
# right_range  = right daughter species' range 
# indep_rate   = density-independent rate for extinction (> 0)
# n_vec_lower  = vector describing the number of of species in each region at time 0
# timewin      = number of timebins

# range: A = c(1), B = c(2), C = c(3), AB = c(1,2), AC = c(1,3), BC = c(2,3)

# test <- heat_map_b(-1,2,0.01,c(3),c(1,2),c(1,5,10),50)

# test <- heat_map_b(0,2,0.01,c(1,2),c(3),c(1,5,10),50) #base case: density independent 

heat_map_b <- function(rate_scale,base_rate,indep_rate,left_range,right_range,
                       n_vec_lower,timewin){
  ##############
  # CONDITION CHECK  #
  ##############
  #
  if(any(n_vec_lower < 0)){
    error("Number of species in any region cannot be negative")
  }
  # browser()
  #
  # Create number of species in each region over time
  n_list <- list()
  for (i in 1:length(n_vec_lower)){
      n_list[[i]]    <- seq(n_vec_lower[i],n_vec_lower[i]+timewin-1,by=1)
  }
  #
  # Create grid for number of species in each range over time
  if (identical(left_range,c(1)) && identical(right_range,c(2))){ #split {A}|{B}
    n_grid <- expand.grid(n_left = n_list[[1]], n_right = n_list[[2]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_A = n_list[[1]],n_B = n_list[[2]])
    n_reg_grid$n_C <- rep(n_list[[3]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(2)) && identical(right_range,c(1))){ #split {B}|{A}
    n_grid <- expand.grid(n_left = n_list[[2]], n_right = n_list[[1]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_B = n_list[[2]],n_A = n_list[[1]])
    n_reg_grid$n_C <- rep(n_list[[3]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(1)) && identical(right_range,c(3))){ #split {A}|{C}
    n_grid <- expand.grid(n_left = n_list[[1]], n_right = n_list[[3]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_A = n_list[[1]],n_C = n_list[[3]])
    n_reg_grid$n_B <- rep(n_list[[2]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(3)) && identical(right_range,c(1))){ #split {C}|{A}
    n_grid <- expand.grid(n_left = n_list[[3]], n_right = n_list[[1]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_C = n_list[[3]],n_A = n_list[[1]])
    n_reg_grid$n_B <- rep(n_list[[2]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(2)) && identical(right_range,c(3))){ #split {B}|{C}
    n_grid <- expand.grid(n_left = n_list[[2]], n_right = n_list[[3]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_B = n_list[[2]],n_C = n_list[[3]])
    n_reg_grid$n_A <- rep(n_list[[1]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(3)) && identical(right_range,c(2))){ #split {C}|{B}
    n_grid <- expand.grid(n_left = n_list[[3]], n_right = n_list[[2]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_C = n_list[[3]],n_B = n_list[[2]])
    n_reg_grid$n_A <- rep(n_list[[1]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(1,2)) && identical(right_range,c(3))){ #split {A,B}|{C}
    n_grid <- expand.grid(n_left = n_list[[1]]+n_list[[2]], n_right = n_list[[3]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_A = n_list[[1]],n_C = n_list[[3]])
    n_reg_grid$n_B <- rep(n_list[[2]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(3)) && identical(right_range,c(1,2))){ #split {C}|{A,B}
    n_grid <- expand.grid(n_left = n_list[[3]], n_right = n_list[[1]]+n_list[[2]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_C = n_list[[3]],n_A = n_list[[1]])
    n_reg_grid$n_B <- expand.grid(n_C = n_list[[3]],n_B = n_list[[2]])$n_B
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(1,3)) && identical(right_range,c(2))){ #split {A,C}|{B}
    n_grid <- expand.grid(n_left = n_list[[1]]+n_list[[3]], n_right = n_list[[2]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_A = n_list[[1]],n_B = n_list[[2]])
    n_reg_grid$n_C <- rep(n_list[[3]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(2)) && identical(right_range,c(1,3))){ #split {B}|{A,C}
    n_grid <- expand.grid(n_left = n_list[[2]], n_right = n_list[[1]]+n_list[[3]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_B = n_list[[2]],n_A = n_list[[1]])
    n_reg_grid$n_C <- expand.grid(n_B = n_list[[2]],n_C = n_list[[3]])$n_C
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(2,3)) && identical(right_range,c(1))){ #split {B,C}|{A}
    n_grid <- expand.grid(n_left = n_list[[2]]+n_list[[3]], n_right = n_list[[1]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_B = n_list[[2]],n_A = n_list[[1]])
    n_reg_grid$n_C <- rep(n_list[[3]],timewin)
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
    #
  } else if (identical(left_range,c(1)) && identical(right_range,c(2,3))){ #split {A}|{B,C}
    n_grid <- expand.grid(n_left = n_list[[1]], n_right = n_list[[2]]+n_list[[3]])
    # Create grid for number of species in each region over time (*should correspond with n_grid)
    n_reg_grid     <- expand.grid(n_A = n_list[[1]],n_B = n_list[[2]])
    n_reg_grid$n_C <- expand.grid(n_A = n_list[[2]],n_C = n_list[[3]])$n_C
    n_reg_grid     <- n_reg_grid[c("n_A","n_B","n_C")]#re-swapping the column order
  } 
 
  # browser()
  # create a list that stores edge weights over time for each model
  edge_weight_log        <- list() # list of edge weights for log model
  edge_weight_loglog     <- list() # list of edge weights for loglog model
  # loop over rows in n_grid
  for(i in 1:length(n_grid$n_left)){
    # initialize empty matrix 
    # store edge weights between two regions
    edge_weight_log[[i]]        <- matrix(NA,nrow = 3,ncol =3)
    edge_weight_loglog[[i]]     <- matrix(NA,nrow = 3,ncol =3)
    # assign the diagonal elements to be equal to 0
    diag(edge_weight_log[[i]])        <- 0
    diag(edge_weight_loglog[[i]])     <- 0
    # compute the edge weight between two regions at each time point (off-diagonal)
    for (j in 1:nrow(edge_weight_log[[i]])){
      for (k in 1:ncol(edge_weight_log[[i]])){
        if(j != k){ # off-diagonal element
          edge_weight_log[[i]][j,k]      <- exp(rate_scale*log(n_reg_grid[[j]][i]+n_reg_grid[[k]][i]-1))
          edge_weight_loglog[[i]][j,k]   <- exp(rate_scale*log(log(n_reg_grid[[j]][i]*n_reg_grid[[k]][i])+1))
        } else if (j == k){ # on-diagonal
          edge_weight_log[[i]][j,k]      <- 0
          edge_weight_loglog[[i]][j,k]   <- 0
        }
      }
    }
  }
  #
  # browser()
  # append density dependent effect and rate over time 
  den_vec_log     <- c()  # log density effect
  den_vec_loglog  <- c()  # loglog density effect
  #
  bkl_vec_log      <- c()  # log density-dependent between-region speciation ranges k,l
  bkl_vec_loglog   <- c()  # loglog density-dependent between-region speciation ranges k,k
  
  # browser()
  
  for (i in 1:length(n_grid$n_left)){
    # compute the range score at each time point (for each model variant)
    range_score_log      <- 0
    range_score_loglog   <- 0
    for (j in 1:length(left_range)){
      for (k in 1:length(right_range)){
          range_score_log      <- range_score_log + 1/((edge_weight_log[[i]][left_range[j],right_range[k]]+edge_weight_log[[i]][right_range[k],left_range[j]])/2)
          range_score_loglog   <- range_score_loglog + 1/((edge_weight_loglog[[i]][left_range[j],right_range[k]]+edge_weight_loglog[[i]][right_range[k],left_range[j]])/2)
      }
    }
    #
    # browser()
    den_vec_log[i]     <- (length(left_range)*length(right_range))*(1/range_score_log) # normalize the density effect for multiplicative model so it return 1 if no density effect
    den_vec_loglog[i]  <- (length(left_range)*length(right_range))*(1/range_score_loglog)
    #
    bkl_vec_log[i]    <- base_rate*indep_rate * den_vec_log[i]
    bkl_vec_loglog[i] <- base_rate*indep_rate * den_vec_loglog[i]
  }
  
  # browser()
  ##############
  # DATA PROCESSING   #
  ##############
  dat_log             <- data.frame(func = rep("log",length(n_grid$n_left)),num_range_left = n_grid[,1], num_range_right = n_grid[,2],
                                    rate_scale = rep(rate_scale,length(n_grid$n_left)),
                                    indep_rate = rep(indep_rate,length(n_grid$n_left)),
                                    density_effect  = den_vec_log,
                                    overall_rate = bkl_vec_log)
  
  
  # browser()
  dat_loglog          <- data.frame(func = rep("loglog",length(n_grid$n_left)),num_range_left = n_grid[,1], num_range_right = n_grid[,2],
                                    rate_scale = rep(rate_scale,length(n_grid$n_left)),
                                    indep_rate = rep(indep_rate,length(n_grid$n_left)),
                                    density_effect  = den_vec_loglog,
                                    overall_rate = bkl_vec_loglog)
  
  #
  # browser()
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog)
  #
  #
  ##############
  # PLOTTINGS  #
  ##############
  func_names  <- c("log","loglog")
  # range_names <- c("A","B","C","AB","AC","BC")
  plot_list = list()
  title_list <- c(paste("log between-region speciation rate\n","rate scale =",rate_scale,"\n",
                        "base rate=",base_rate,",","independent rate=",indep_rate),
                  paste("log-log between-region speciation rate\n","rate scale =",rate_scale,"\n",
                        "base rate=",base_rate,",","independent rate=",indep_rate))
  # make plot for each model
  # browser()
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    # browser()
    # browser()
    #
    p = ggplot(dat_i, aes(x=num_range_left ,y=num_range_right,fill=overall_rate))
    #
    p = p+geom_tile()
    # plot labels 
    p = p + xlab("#species in left range") + 
        ylab("#species in right range") + 
        ggtitle(title_list[i]) + theme(plot.title = element_text(hjust = 0.5))
    #
    plot_list[[i]] = p
    print(p)
  }
  plot_list[[1]]
  plot_list[[2]]
  #
  # return(dat_comb)
  return(n_reg_grid)
}