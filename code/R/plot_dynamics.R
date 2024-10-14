library(ggplot2)

##############
# EXTINCTION  #
##############

# t_lower     = lower bound of runtime t
# t_upper     = upper bound of runtime t
# e_d         = extinction rate factor in home location (-inf,inf)
# base        = base ext rate
# base_lower  = lower bound of density-independent rate (> 0 )
# base_upper  = upper bound of density-independent rate (> 0)
# n_lower.    = minimum number of species in i
# type        = species count dynamics 
#               ("constant increase","constant decrease",
#               "random","constant")
# timewin     = number of timebins
# opt         = turn off (0) or on (1) "(ni,nj)" text on plots

## Reduce to base case (density-independent rate) with base rate = 0.1

# test <- plot_e(0,5,0,0.1,0.1,10,"constant increase",100,0) # set w_d = 0 (i.e., turning off density feature)
# test <- plot_e(0,5,0.5,0.1,0.1,1,"constant",100,0) # set n_i(t) = 1 (i.e., only have one species)

# test <- plot_e(0,500,0.05,0.1,0.1,10,"constant increase",10000,0) # density effect increases extinction rate
# test <- plot_e(0,500,-0.05,0.1,0.1,10,"constant increase",10000,0) # density effect decreases extinction rate

plot_e <- function(t_lower,t_upper,e_d,base,base_lower,base_upper,n_lower,type,timewin,opt){
  ##############
  # CONDITION CHECK  #
  ##############
  
  # check if input if satisfied
  if(base_lower <= 0){
    error("Density independent rate must be postive")
  }
  #
  if(n_lower < 0){
    error("Number of species cannot be negative")
  }
  #
  if(t_lower < 0){
    error("Minimum runtime cannot be negative")
  }
  #
  ##############
  # SETTINGS   #
  ##############
  
  # density independent extinction rate
  e_i_i   <- seq(base_lower,base_upper,length.out = timewin)
  # extinction rate factor 
  e_d     <- e_d
  # time span
  tspan   <- seq(t_lower,t_upper,length.out = timewin)
  # number of species over time
  if (type == "constant increase"){
    n_i     <- seq(n_lower,n_lower+timewin-1,by=1)
  } else if (type == "constant decrease"){
    n_i     <- sort(seq(n_lower,n_lower+timewin-1,by=1),decreasing = TRUE)
  } else if (type == "random"){
    n_i     <- sample.int(100, timewin)
  } else if (type == "constant"){
    n_i     <- rep(n_lower,timewin)
  }
  # append density dependent effect and rate over time 
  # den_vec_lin  <- c()  # additive density effect
  den_vec_log     <- c()  # log density effect
  den_vec_loglog  <- c()  # loglog density effect
  #
  # ei_vec_lin   <- c()  # linear density-dependent extinction
  ei_vec_log     <- c()  # log density-dependent extinction
  ei_vec_loglog  <- c()  # loglog density-dependent extinction
  
  for (i in 1:length(tspan)){
    # den_vec_lin[i]  <- (n_i[i]-1)*e_d
    # den_vec_mult[i] <- exp(e_d)^(n_i[i]-1)
    den_vec_log[i]     <- n_i[i]^e_d # Equivalent to exp(e_d*log(n_i)), linear trend if e_d=1
    # den_vec_log[i]  <- 1 + log(n_i[i])*e_d
    den_vec_loglog[i]  <- exp(e_d*log(log(n_i[i])+1))
    #
    # ei_vec_lin[i]  <- e_i_i[i] + den_vec_lin[i]
    ei_vec_log[i]     <- base*e_i_i[i] * den_vec_log[i]
    ei_vec_loglog[i]  <- base*e_i_i[i] * den_vec_loglog[i]
  }
  # browser()
  
  ##############
  # DATA PROCESSING   #
  ##############
  
  # store data
  # dat_linear        <- data.frame(func = rep("additive",timewin),time = round(tspan,2), num_spec = n_i, rate_factor = rep(e_d,timewin),
  #                                 density_indep_rate = round(e_i_i,3),density_effect = round(den_vec_lin,3), 
  #                                 density_dep_rate = round(ei_vec_lin,3))
  #
  dat_log             <- data.frame(func = rep("log",timewin),time = round(tspan,2), num_spec = n_i, rate_factor = rep(e_d,timewin),
                                  density_indep_rate = round(e_i_i,3),density_effect = round(den_vec_log,3), 
                                  density_dep_rate = round(ei_vec_log,3))
  #
  dat_loglog          <- data.frame(func = rep("loglog",timewin),time = round(tspan,2), num_spec = n_i, rate_factor = rep(e_d,timewin),
                                  density_indep_rate = round(e_i_i,3),density_effect = round(den_vec_loglog,3), 
                                  density_dep_rate = round(ei_vec_loglog,3))
  
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog,by="time") 
  dat_comb <- dat_comb[-(2*timewin+1),]
#
  dat_comb$time <- as.numeric(dat_comb$time)
  dat_comb$density_dep_rate <- as.numeric(dat_comb$density_dep_rate)
# order by time  
  dat_comb <- dat_comb[order(dat_comb$time),]
  rownames(dat_comb) <- NULL
  #
  # dat_filter <- subset(dat,dat$density_dep_rate >= 0)
  # browser()
  ##############
  # PLOTTINGS  #
  ##############
  func_names <- c("log","loglog")
  plot_list = list()
  title_list <- c("standard extinction","logarithmic extinction")
  # make plot
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    rownames(dat_i) <- NULL
    p = ggplot(dat_i, aes(x=time, y=density_dep_rate) )
    #
    p = p+geom_point()
    #
    p = p+geom_line()
    #
    p = p + geom_hline(yintercept=0, linetype="dashed", color = "red")
    #
    if (opt == 1){
      p = p + geom_text(aes(label = num_spec),
                        hjust = -0.5, size = 2)
    } else if (opt == 0){
      p = p
    }
    # plot labels 
    p = p + xlab("time") + ylab("density dependent extinction rate") + ggtitle(title_list[i])
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
# WITHIN-REGION SPECIATION  #
##############

# t_lower     = lower bound of runtime t
# t_upper     = upper bound of runtime t
# w_d         = within-region rate factor in home location (-inf,inf)
# base        = base w rate
# base_lower  = lower bound of density-independent rate (> 0 )
# base_upper  = upper bound of density-independent rate (> 0)
# n_lower.    = minimum number of species in i
# type        = species count dynamics 
#               ("constant increase","constant decrease",
#               "random","constant")
# timewin     = number of timebins
# opt         = turn off (0) or on (1) "(ni,nj)" text on plots

## Reduce to base case (density-independent rate) with base rate = 0.1
# test <- plot_w(0,5,0,0.1,0.1,10,"constant increase",100,0) # set w_d = 0 (i.e., turning off density feature)
# test <- plot_w(0,5,0.5,0.1,0.1,1,"constant",100,0) # set n_i(t) = 1 (i.e., only have one species)

# test <- plot_w(0,500,0.05,0.1,0.1,0,"constant increase",10000,1) #density effect increases within-region speciation rate 

# test <- plot_w(0,500,-0.05,0.1,0.1,0,"constant increase",10000,1)  #density effect decreases within-region speciation rate 

plot_w <- function(t_lower,t_upper,w_d,base,base_lower,base_upper,n_lower,type,timewin,opt){
  ##############
  # CONDITION CHECK  #
  ##############
  
  # check if input if satisfied
  if(base_lower <= 0){
    error("Density independent rate must be postive")
  }
  #
  if(n_lower < 0){
    error("Number of species cannot be negative")
  }
  #
  if(t_lower < 0){
    error("Minimum runtime cannot be negative")
  }
  #
  ##############
  # SETTINGS   #
  ##############
  
  # density independent extinction rate
  w_i_i   <- seq(base_lower,base_upper,length.out = timewin)
  # extinction rate factor 
  w_d     <- w_d
  # time span
  tspan   <- seq(t_lower,t_upper,length.out = timewin)
  # number of species over time
  if (type == "constant increase"){
    n_i     <- seq(n_lower,n_lower+timewin-1,by=1)
  } else if (type == "constant decrease"){
    n_i     <- sort(seq(n_lower,n_lower+timewin-1,by=1),decreasing = TRUE)
  } else if (type == "random"){
    n_i     <- sample.int(100, timewin)
  } else if (type == "constant"){
    n_i     <- rep(n_lower,timewin)
  }
  # append density dependent effect and rate over time 
  # den_vec_lin  <- c()  # additive density effect
  den_vec_log     <- c()  # log density effect
  den_vec_loglog  <- c()  # log-log density effect
  #
  # wi_vec_lin   <- c()  # linear density-dependent within-region speciation
  wi_vec_log      <- c()  # log density-dependent within-region speciation
  wi_vec_loglog   <- c()  # log-log density-dependent within-region speciation
  
  for (i in 1:length(tspan)){
    # den_vec_lin[i]  <- (n_i[i]-1)*w_d
    # den_vec_mult[i] <- exp(w_d)^(n_i[i]-1)
    den_vec_log[i]     <- exp(w_d*log(n_i[i]))
    # den_vec_log[i]  <- 1+log(n_i[i])*w_d
    den_vec_loglog[i]  <- exp(w_d*log(log(n_i[i])+1))
    #
    # wi_vec_lin[i]  <- w_i_i[i] + den_vec_lin[i]
    wi_vec_log[i]     <- base*w_i_i[i] * den_vec_log[i]
    wi_vec_loglog[i]  <- base*w_i_i[i] * den_vec_loglog[i]
  }
  # browser()
  
  ##############
  # DATA PROCESSING   #
  ##############
  
  # store data
  # dat_linear        <- data.frame(func = rep("additive",timewin),time = round(tspan,2), num_spec = n_i, rate_factor = rep(w_d,timewin),
  #                                 density_indep_rate = round(w_i_i,3),density_effect = round(den_vec_lin,3), 
  #                                 density_dep_rate = round(wi_vec_lin,3))
  #
  dat_log              <- data.frame(func = rep("log",timewin),time = round(tspan,2), num_spec = n_i, rate_factor = rep(w_d,timewin),
                                  density_indep_rate = round(w_i_i,3),density_effect = round(den_vec_log,3), 
                                  density_dep_rate = round(wi_vec_log,3))
  #
  dat_loglog           <- data.frame(func = rep("loglog",timewin),time = round(tspan,2), num_spec = n_i, rate_factor = rep(w_d,timewin),
                                  density_indep_rate = round(w_i_i,3),density_effect = round(den_vec_loglog,3), 
                                  density_dep_rate = round(wi_vec_loglog,3))
  # browser()
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog,by="time") 
  dat_comb <- dat_comb[-(2*timewin+1),]
  #
  dat_comb$time <- as.numeric(dat_comb$time)
  dat_comb$density_dep_rate <- as.numeric(dat_comb$density_dep_rate)
  # order by time  
  dat_comb <- dat_comb[order(dat_comb$time),]
  rownames(dat_comb) <- NULL
  #
  # dat_filter <- subset(dat,dat$density_dep_rate >= 0)
  # browser()
  ##############
  # PLOTTINGS  #
  ##############
  func_names <- c("log","loglog")
  plot_list = list()
  title_list <- c("standard within-region speciation","logarithmic within-region speciation")
  # make plot
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    rownames(dat_i) <- NULL
    p = ggplot(dat_i, aes(x=time, y=density_dep_rate) )
    #
    p = p+geom_point()
    #
    p = p+geom_line()
    #
    p = p + geom_hline(yintercept=0, linetype="dashed", color = "red")
    #
    if (opt == 1){
      p = p + geom_text(aes(label = num_spec),
                        hjust = -0.5, size = 2)
    } else if (opt == 0){
      p = p
    }
    #
    # plot labels 
    p = p + xlab("time") + ylab("density dependent within-region speciation rate") + ggtitle(title_list[i])
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
# DISPERSAL  #
##############

# t_lower     = lower bound of runtime t
# t_upper     = upper bound of runtime t
# d_d_i       = dispersal rate factor in home location (-inf,inf)
# d_d_j       = dispersal rate factor in target location (inf,inf)
# base        = base dispersal rate
# base_lower  = lower bound of density-independent rate (> 0 )
# base_upper  = upper bound of density-independent rate (> 0)
# n_i_lower.  = minimum number of species in home location (> 0)
# n_j_lower.  = minimum number of species in target location (>= 0)
# type        = species count dynamics 
                # ("constant increase both","constant decrease both","constant decrease home increase away","constant increase home decrease away",
                #  "random both","constant")
# timewin     = number of timebins
# opt         = turn off (0) or on (1) "(ni,nj)" text on plots

## Reduce to base case (density-independent rate) with base rate = 0.1
# test <- plot_d_i_j(0,5,0,0,0.1,0.1,10,10,"constant increase both",100,1) #set d_d_i = d_d_j = 0 (i.e. turning off density feature)
# test <- plot_d_i_j(0,5,0.05,0.05,0.1,0.1,1,0,"constant",100,1) #set n_i(t) = 1 and n_j(t) = 0 (i.e. only 1 species at home, and no species at target)
# test <- plot_d_i_j(0,5,0,1,0.1,0.1,10,0,"constant",100,1) #set d_d_i = 0 and n_j(t) = 0 (i.e. turning off density feature at home, and no species in j)
# test <- plot_d_i_j(0,5,10,0,0.1,0.1,1,10,"constant",100,1) #set n_i(t) = 1 and d_d_j = 0 (i.e. only 1 species at home, and turning off density feature in target)

# test <- plot_d_i_j(0,500,0.05,0.05,0.1,0.1,10,10,"constant increase both",10000,0) # density effect increases dispersal rate i to j with increasing population 
                                                                                # (i.e. d_i_i & d_i_j > 0)

# test <- plot_d_i_j(0,500,0.01,-0.05,0.1,0.1,10,10,"constant increase both",10000,0) # DD in j overpowers DD in i to decrease the rate
# test <- plot_d_i_j(0,500,-0.01,0.05,0.1,0.1,10,10,"constant increase both",10000,0) # DD in j overpowers DD in i to increase the rate

# test <- plot_d_i_j(0,500,0.05,-0.01,0.1,0.1,10,10,"constant increase both",10000,0) # DD in i overpowers DD in j to increase the rate
# test <- plot_d_i_j(0,500,-0.05,0.01,0.1,0.1,10,10,"constant increase both",10000,0) # DD in i overpowers DD in j to decrease the rate

plot_d_i_j <- function(t_lower,t_upper,d_d_i,d_d_j,base,base_lower,base_upper,n_i_lower,n_j_lower,type,timewin,opt){
  ##############
  # CONDITION CHECK  #
  ##############
  
  # check if input if satisfied
  if(base_lower <= 0){
    error("Density independent rate must be postive")
  }
  #
  if(n_i_lower <= 0){
    error("Number of species in home region must be positive")
  }
  #
  if(n_j_lower < 0){
    error("Number of species in target region must be non-negative")
  }
  #
  if(t_lower < 0){
    error("Minimum runtime cannot be negative")
  }
  #
  ##############
  # SETTINGS   #
  ##############
  
  # density independent dispersal rate i to j 
  d_i_j        <- seq(base_lower,base_upper,length.out = timewin)
  # dispersal rate factor in home location i 
  d_d_home     <- d_d_i
  # dispersal rate factor in target location j
  d_d_target   <- d_d_j
  # time span
  tspan   <- seq(t_lower,t_upper,length.out = timewin)
  # number of species over time in i and j
  if (type == "constant increase both"){
    n_i     <- seq(n_i_lower,n_i_lower+timewin-1,by=1)
    n_j     <- seq(n_j_lower,n_j_lower+timewin-1,by=1)
  } else if (type == "constant decrease both"){
    n_i     <- sort(seq(n_i_lower,n_i_lower+timewin-1,by=1),decreasing = TRUE)
    n_j     <- sort(seq(n_j_lower,n_j_lower+timewin-1,by=1),decreasing = TRUE)
  } else if (type == "constant decrease home increase away"){
    n_i     <- sort(seq(n_i_lower,n_i_lower+timewin-1,by=1),decreasing = TRUE)
    n_j     <- sort(seq(n_j_lower,n_j_lower+timewin-1,by=1),decreasing = FALSE)
  } else if (type == "constant increase home decrease away"){
    n_i     <- sort(seq(n_i_lower,n_i_lower+timewin-1,by=1),decreasing = FALSE)
    n_j     <- sort(seq(n_j_lower,n_j_lower+timewin-1,by=1),decreasing = TRUE)
  } else if (type == "random both"){
    n_i     <- sample.int(100, timewin)
    n_j     <- sample.int(100, timewin)
  } else if (type == "constant"){
    n_i     <- rep(n_i_lower,timewin)
    n_j     <- rep(n_j_lower,timewin)
  }
  # browser()
  # append density dependent effect and rate over time 
  # den_vec_lin    <- c()  # additive density effect
  den_vec_log       <- c()  # log density effect
  den_vec_loglog    <- c()  # loglog density effect
  #
  # dij_vec_lin     <- c()  # linear density-dependent dispersal rate i to j
  dij_vec_log        <- c()  # log density-dependent dispersal rate i to j
  dij_vec_loglog     <- c()  # loglog density-dependent dispersal rate i to j
  
  for (i in 1:length(tspan)){
    # den_vec_lin[i]     <- (n_i[i]-1)*d_d_home + (n_j[i])*d_d_target
    # den_vec_mult[i]  <- (1+log(n_i[i])*d_d_home+log(n_j[i]+1)*d_d_target)*log(exp(1)+(d_d_home*d_d_target*n_i[i]*n_j[i])/(1+d_d_home*d_d_target))
    # den_vec_mult[i]    <- exp(d_d_home)^(n_i[i]-1)*exp(d_d_target)^(n_j[i])
    den_vec_log[i]     <- exp(d_d_i*log(n_i[i])+d_d_j*log(n_j[i]+1))
    den_vec_loglog[i]  <- exp(d_d_i*log(log(n_i[i])+1)+d_d_j*log(log(n_j[i]+1)+1))
    #
    # dij_vec_lin[i]    <- d_i_j[i] + den_vec_lin[i]
    dij_vec_log[i]       <- base*d_i_j[i] * den_vec_log[i]
    dij_vec_loglog[i]    <- base*d_i_j[i] * den_vec_loglog[i]
  }
  # browser()
  
  ##############
  # DATA PROCESSING   #
  ##############
  
  # store data
  # dat_linear        <- data.frame(func = rep("additive",timewin),time = round(tspan,2), num_spec_i = n_i,num_spec_j = n_j,
  #                                 rate_factor_i = rep(d_d_home,timewin),rate_factor_j = rep(d_d_target,timewin),
  #                                 density_indep_rate = round(d_i_j,3),density_effect = round(den_vec_lin,3), 
  #                                 density_dep_rate = round(dij_vec_lin,3))
  #
  dat_log             <- data.frame(func = rep("log",timewin),time = round(tspan,2), num_spec_i = n_i,num_spec_j = n_j,
                                    rate_factor_i =rep(d_d_home,timewin),rate_factor_j = rep(d_d_target,timewin),
                                    density_indep_rate = round(d_i_j,3),density_effect = round(den_vec_log,3), 
                                    density_dep_rate = round(dij_vec_log,3))
  #
  dat_loglog          <- data.frame(func = rep("loglog",timewin),time = round(tspan,2), num_spec_i = n_i,num_spec_j = n_j,
                                  rate_factor_i =rep(d_d_home,timewin),rate_factor_j = rep(d_d_target,timewin),
                                  density_indep_rate = round(d_i_j,3),density_effect = round(den_vec_loglog,3), 
                                  density_dep_rate = round(dij_vec_loglog,3))
  #
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog,by="time") 
  dat_comb <- dat_comb[-(2*timewin+1),]
  #
  dat_comb$time <- as.numeric(dat_comb$time)
  dat_comb$density_dep_rate <- as.numeric(dat_comb$density_dep_rate)
  # order by time  
  dat_comb <- dat_comb[order(dat_comb$time),]
  rownames(dat_comb) <- NULL
  #
  # dat_filter <- subset(dat,dat$density_dep_rate >= 0)
  # browser()
  ##############
  # PLOTTINGS  #
  ##############
  func_names <- c("log","loglog")
  plot_list = list()
  title_list <- c("standard dispersal","logarithmic dispersal")
  # make plot
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    rownames(dat_i) <- NULL
    p = ggplot(dat_i, aes(x=time, y=density_dep_rate) )
    #
    p = p+geom_point()
    #
    p = p+geom_line()
    #
    p = p + geom_hline(yintercept=0, linetype="dashed", color = "red")
    #
    if (opt == 1){ # plot the text 
      p = p + geom_text(aes(label = paste("(",num_spec_i,",",num_spec_j,")")),
                        hjust = -0.5, size = 2)
    } else if (opt ==0){ # dont plot the text
      p = p
    }
    #
    # plot labels 
    p = p + xlab("time") + ylab("density dependent dispersal rate i to j") + ggtitle(title_list[i])
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
# BETWEEN-REGION SPECIATION  #
####
##############

# t_lower         = lower bound of runtime t
# t_upper         = upper bound of runtime t

# left range      = range of left species, define as set of regions where species occupy, = c(region,region,...)
# right range     = range of right species, define as set of regions where species occupy, = c(region,region,...)
#* must be in order: i.e. left range = c(1,2) is valid, but c(2,1) is not
# 1 = region A, 2 = region B, 3 = region C

# b_d             = between-region rate factor
# base            = base B rate
# base_lower      = lower bound of density-independent rate
# base_upper      = upper bound of density-independent rate

# Setting base_lower == base_upper means we have constant density-independent rate 

# n_lower_vec     = vector of minimum number of species in each region

# type            = species count dynamics ("constant increase","constant decrease","random","both")
# timewin         = number of timebins
# opt             = turn off (0) or on (1) "(ni,nj)" text on plots

## Reduce to base case (density-independent rate) with base rate = 0.1
# test <- plot_b(0,5,c(1),c(3),0,0.1,0.1,c(1,1,1),"constant decrease",100,0)

# test <- plot_b(0,500,c(1),c(3),0.05,0.1,0.1,c(1,1,1),"constant increase",10000,0) # density effect increases between region speciation ranges k and l as species increases
# test <- plot_b(0,500,c(1),c(3),-0.05,0.1,0.1,c(1,1,1),"constant increase",10000,0) # density effect decreases between region speciation ranges k and l as species increases

plot_b <- function(t_lower,t_upper,left_range,right_range,b_d,base,base_lower,base_upper,n_vec_lower,type,timewin,opt){
  ##############
  # CONDITION CHECK  #
  ##############
  # check if input if satisfied
  if(base_lower <= 0){
    error("Density independent rate must be postive")
  }
  #
  if(any(n_vec_lower < 0)){
    error("Number of species in any region cannot be negative")
  }
  #
  if(t_lower < 0){
    error("Minimum runtime cannot be negative")
  }
  #
  ##############
  # SETTINGS   #
  ##############
  # browser()
  # density independent extinction rate
  b_k_l   <- seq(base_lower,base_upper,length.out = timewin)
  # between-region rate factor 
  b_d     <- b_d
  # time span
  tspan   <- seq(t_lower,t_upper,length.out = timewin)
  # number of species over time in each region 
  # browser()
  n_list <- list()
  for (i in 1:length(n_vec_lower)){
    if (type == "constant increase"){
      n_list[[i]]    <- seq(n_vec_lower[i],n_vec_lower[i]+timewin-1,by=1)
    } else if (type == "constant decrease"){
      n_list[[i]]    <- sort(seq(n_vec_lower[i],n_vec_lower[i]+timewin-1,by=1),decreasing = TRUE)
    } else if (type == "random"){
      n_list[[i]]    <- sample.int(100, timewin)
    } else if (type == "constant"){
      n_list[[i]]     <- rep(n_vec_lower[i],timewin)
    }
  }
  # create a list that stores edge weights over time for each model
  # edge_weight_add     <- list() # list of edge weights for additive model
  edge_weight_log        <- list() # list of edge weights for log model
  edge_weight_loglog     <- list() # list of edge weights for loglog model
  # looping for every time bin
  for(i in 1:length(tspan)){
    # initialize empty matrix at t_i
    # edge_weight_add[[i]]     <- matrix(NA,nrow = 3,ncol =3)
    edge_weight_log[[i]]        <- matrix(NA,nrow = 3,ncol =3)
    edge_weight_loglog[[i]]     <- matrix(NA,nrow = 3,ncol =3)
    # assign the diagonal elements to be equal to 0
    # diag(edge_weight_add[[i]])     <- 0
    diag(edge_weight_log[[i]])   <- 0
    diag(edge_weight_loglog[[i]])     <- 0
    # compute the edge weight between two regions at each time point
    for (j in 1:nrow(edge_weight_log[[i]])){
      for (k in 1:ncol(edge_weight_log[[i]])){
        if (j != k){
          edge_weight_log[[i]][j,k]      <- exp(b_d*log(n_list[[j]][i]+n_list[[k]][i]-1))
          edge_weight_loglog[[i]][j,k]   <- exp(b_d*log(log(n_list[[j]][i]*n_list[[k]][i])+1))
        } else if (j == k){
          edge_weight_log[[i]][j,k]      <- 0
          edge_weight_loglog[[i]][j,k]   <- 0
        }
      }
    }
  }
  # browser()
  #
  # browser()
  # append density dependent effect and rate over time 
  # den_vec_lin  <- c()  # additive density effect
  den_vec_log     <- c()  # log density effect
  den_vec_loglog  <- c()  # loglog density effect
  #
  # bkl_vec_lin   <- c()  # linear density-dependent between-region speciation ranges k,l
  bkl_vec_log      <- c()  # log density-dependent between-region speciation ranges k,l
  bkl_vec_loglog   <- c()  # loglog density-dependent between-region speciation ranges k,k

  for (i in 1:length(tspan)){
    # compute the range score at each time point (for each model variant)
    # range_score_add   <- 0
    range_score_log <- 0
    range_score_loglog   <- 0
    for (j in 1:length(left_range)){
      for (k in 1:length(right_range)){
        range_score_log      <- range_score_log + 1/((edge_weight_log[[i]][left_range[j],right_range[k]]+edge_weight_log[[i]][right_range[k],left_range[j]])/2)
        range_score_loglog   <- range_score_loglog + 1/((edge_weight_loglog[[i]][left_range[j],right_range[k]]+edge_weight_loglog[[i]][right_range[k],left_range[j]])/2)
      }
    }
    #
    # browser()
    # den_vec_lin[i]  <- b_d * 1/range_score_add
    den_vec_log[i]     <- (length(left_range)*length(right_range))*(1/range_score_log) # normalize the density effect for multiplicative model so it return 1 if no density effect
    den_vec_loglog[i]  <- (length(left_range)*length(right_range))*(1/range_score_loglog)
    #
    # bkl_vec_lin[i]  <- b_k_l[i] + den_vec_lin[i]
    bkl_vec_log[i]    <- base*b_k_l[i] * den_vec_log[i]
    bkl_vec_loglog[i] <- base*b_k_l[i] * den_vec_loglog[i]
  }
  # browser()
  
  ##############
  # DATA PROCESSING   #
  ##############
  
  # store data
  # dat_linear        <- data.frame(func = rep("additive",timewin),time = round(tspan,2), num_spec_A =  n_list[[1]], 
  #                                 num_spec_B =  n_list[[2]],  num_spec_C =  n_list[[3]],rate_factor = rep(b_d,timewin),
  #                                 density_indep_rate = round(b_k_l,3),density_effect = round(den_vec_lin,3), 
  #                                 density_dep_rate = round(bkl_vec_lin,3))
  #
  dat_log              <- data.frame(func = rep("log",timewin),time = round(tspan,2), num_spec_A =  n_list[[1]], 
                                  num_spec_B =  n_list[[2]],  num_spec_C =  n_list[[3]], rate_factor = rep(b_d,timewin),
                                  density_indep_rate = round(b_k_l,3),density_effect = round(den_vec_log,3), 
                                  density_dep_rate = round(bkl_vec_log,3))
  #
  dat_loglog           <- data.frame(func = rep("loglog",timewin),time = round(tspan,2), num_spec_A =  n_list[[1]], 
                                  num_spec_B =  n_list[[2]],  num_spec_C =  n_list[[3]], rate_factor = rep(b_d,timewin),
                                  density_indep_rate = round(b_k_l,3),density_effect = round(den_vec_loglog,3), 
                                  density_dep_rate = round(bkl_vec_loglog,3))
  # combined list (log and loglog)
  dat_comb <- rbind(dat_log,dat_loglog,by="time") 
  dat_comb <- dat_comb[-(2*timewin+1),]
  #
  dat_comb$time <- as.numeric(dat_comb$time)
  dat_comb$density_dep_rate <- as.numeric(dat_comb$density_dep_rate)
  # order by time  
  dat_comb <- dat_comb[order(dat_comb$time),]
  rownames(dat_comb) <- NULL
  #
  # dat_filter <- subset(dat,dat$density_dep_rate >= 0)
  # browser()
  ##############
  # PLOTTINGS  #
  ##############
  func_names <- c("log","loglog")
  plot_list = list()
  title_list <- c("standard between-region","logarithmic between-region")
  # make plot
  for (i in 1:2){
    dat_i <- dat_comb[dat_comb$func==func_names[i],]
    rownames(dat_i) <- NULL
    p = ggplot(dat_i, aes(x=time, y=density_dep_rate) )
    #
    p = p+geom_point()
    #
    p = p+geom_line()
    #
    p = p + geom_hline(yintercept=0, linetype="dashed", color = "red")
    #
    if (opt == 1){
      p = p + geom_text(aes(label = paste("(",num_spec_A,",",num_spec_B,",",num_spec_C,")")),
                        hjust = -0.5, size = 2)
    } else if (opt == 0){
      p = p
    }
    #
    # plot labels 
    p = p + xlab("time") + ylab("density dependent between-region speciation rate") + ggtitle(title_list[i])
    #
    plot_list[[i]] = p
    print(p)
  }
  plot_list[[1]]
  plot_list[[2]]
  #
  return(dat_comb)
}
