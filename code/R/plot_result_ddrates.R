# Show trajectory under the case of no other factors except density-dependent, and each region share
# equal local equiilibrium capacity (indendent of which region)
#
# type: "simulate" if it is true generating values "estimate" if it is estimated from phyddle analyses
# file_name: name of file being read (only for "individual" summary_type)
# max_species: upper bound on number of species occupying a region
# summary_type: "batch" (take average across all trees and plot) or "individual" (over an individual tree)

#e.g.
# plot_w_e_div("estimate","sim.empirical_est.labels_num.csv",100,1000,summary_type = "individual")

plot_w_e_div <- function(type,file_name = NULL,min_species,max_species,summary_type){
  ##############
  # LOAD R LIBRARIES #
  ##############
  library(gridExtra)
  library(ggplot2)
  #
  ##############
  # FILE SYSTEM #
  ##############
  # browser()
  base_path      = Sys.getenv("HOME") # get to user's base directory
  fp             = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","ddwithin_jl_final_newset_2")
  if (file.exists(fp) == TRUE){
    in_fp          = paste0(fp,'/',type,'/')
  } else {
    warning(paste0("path ",fp, " does not exist"))
  }
  #
  ##############
  # CONDITION CHECK  #
  ##############
  ####################################
  #### Individual vs batch data file ####
  ####################################
  #### Individual dataset 
  if (summary_type == "individual"){
      ####################################
      #### Estimate vs Simulate data file ####
      ####################################
      ##### Based on estimated param values #####
    if (type == "estimate"){
        setwd(in_fp)
        # load the param file 
        param_dat = read.csv(file_name)
        if (file.exists(file_name) == TRUE){
          # read and store the values
          #### Get rate scalar for w
          if ("exp_scalar_w_value" %in% names(param_dat) == TRUE){
            w_d = log(param_dat[1,"exp_scalar_w_value"]) # rate scalar for within speciation
          } else {
            w_d = 0 # must be the case where no density effect on w process 
          }
          #### Get rate scalar for e
          # browser()
          if ("exp_scalar_e_value" %in% names(param_dat) == TRUE){
            e_d = log(param_dat[1,"exp_scalar_e_value"]) # rate scalar for extinction
          } else {
            e_d = 0 # must be the case where no density effect on e process 
          }
          #### Get rate scalar for base rho_w
          if ("log_base_w_value" %in% names(param_dat) == TRUE){
            rho_w = exp(param_dat[1,"log_base_w_value"]) # base rate for within speciation
          } else {
            warning("Invalid data, does not contain within-region speciation process")
          }
          #### Get rate scalar for base rho_e
          if ("log_base_e_value" %in% names(param_dat) == TRUE){
            rho_e = exp(param_dat[1,"log_base_e_value"]) # base rate for extinction
          } else {
            warning("Invalid data, does not contain extinction process")
          }
        } else {
          warning("the file name ", file_name," does not exist in the current estimate directory")
        }
        ##### Based on simulate param values #####
    } else if (type == "simulate"){
        setwd(in_fp)
        # load the param file 
        param_dat = read.csv(file_name)
        if (file.exists(file_name) == TRUE){
          # read and store the values
          #### Get rate scalar for w
          if ("exp_scalar_w" %in% names(param_dat) == TRUE){
            w_d = log(param_dat[1,"exp_scalar_w"]) # rate scalar for within speciation
          } else {
            warning("invalid data, does not contain rate scalar for within speciation")
          }
          #### Get rate scalar for e
          if ("exp_scalar_e" %in% names(param_dat) == TRUE){
            e_d = log(param_dat[1,"exp_scalar_e"]) # rate scalar for extinction
          } else {
            warning("invalid data, does not contain rate scalar for extinction")
          }
          #### Get rate scalar for base rho_w
          if ("log_base_w" %in% names(param_dat) == TRUE){
            rho_w = exp(param_dat[1,"log_base_w"]) # base rate for within speciation
          } else {
            warning("Invalid data, does not contain within-region speciation process")
          }
          #### Get rate scalar for base rho_e
          if ("log_base_e" %in% names(param_dat) == TRUE){
            rho_e = exp(param_dat[1,"log_base_e"]) # base rate for extinction
          } else {
            warning("Invalid data, does not contain extinction process")
          }
        } else {
          warning("the file name ", file_name," does not exist in the current simulate directory")
        }
      } else {
        warning("invalid type, can only plot using true generating parameters or estimated values")
      }
  #### Batch dataset if using simulated data 
  } else if (summary_type == "batch"){
      if (is.null(file_name) == TRUE || is.null(file_name == FALSE)){ # Does not need file_name when summary_type = batch
        if (type == "simulate"){
          setwd(in_fp)
          temp = list.files(pattern = "sim\\.[0-9]+\\.labels\\.csv", full.names = FALSE)
          indices_vec = as.numeric(gsub("sim\\.([0-9]+)\\.labels\\.csv", "\\1", temp)) # get the indices from these files 
          # Initialize vectors to store 
          w_d_vec     = c() # vector of w_d across all simulations
          e_d_vec     = c() # vector of e_d across all simulations
          rho_w_vec   = c() # vector of rho_w across all simulations
          rho_e_vec   = c() # vector of rho_e across all simulations 
          # loop over these files
          for (i in length(indices_vec)){
            # read param data for sim i
            param_dat = read.csv(paste0("sim.",indices_vec[i],".labels.csv"))
            #
            # Get rate scalar for w
            if ("exp_scalar_w_value" %in% names(param_dat) == TRUE){
              w_d = log(param_dat[1,"exp_scalar_w_value"]) # rate scalar for within speciation
            } else {
              w_d = 0 # must be the case where no density effect on w process 
            }
            w_d_vec = append(w_d_vec,w_d) # append w_d from sim_i to w_d_vec
            # Get rate scalar for e
            if ("exp_scalar_e_value" %in% names(param_dat) == TRUE){
              e_d = log(param_dat[1,"exp_scalar_e_value"]) # rate scalar for extinction
            } else {
              e_d = 0 # must be the case where no density effect on e process 
            }
            e_d_vec = append(e_d_vec,e_d) # append e_d from sim_i to e_d_vec
            # Get rate scalar for base rho_w
            if ("log_base_w_value" %in% names(param_dat) == TRUE){
              rho_w = exp(param_dat[1,"log_base_w_value"]) # base rate for within speciation
            } else {
              warning("Invalid data, does not contain within-region speciation process")
            }
            rho_w_vec = append(rho_w_vec,rho_w)
            # Get rate scalar for base rho_e
            if ("log_base_e_value" %in% names(param_dat) == FALSE){
              rho_e = exp(param_dat[1,"log_base_e_value"]) # base rate for extinction
            } else {
              warning("Invalid data, does not contain extinction process")
            }
            rho_e_vec = append(rho_e_vec,rho_e)
          }
          # get the average
          w_d   = mean(w_d_vec)
          e_d   = mean(e_d_vec)
          rho_w = mean(rho_w_vec)
          rho_e = mean(rho_e_vec)
        } else {
          warning("invalid type, can only plot using true generating parameters or estimated values")
        }
      }
  } else {
    warning("invalid summary type, can only take individual or batch datasets")
  }

  ##############
  # Density-dependent overall rate#
  ##############
  # number of species over time
  if (min_species <= 0){
    warning("Minimum number of species cannot be 0 or negative")
  } else {
    n_i = seq(min_species,max_species)
  }
  # Initialize overall speciation rate vector
  spec_rate_log       = c() # Log DDG model
  spec_rate_standard  = c() # Standard DDG model
  # Initialize overall extinction rate vector
  ext_rate_log        = c() # Log DDG model
  ext_rate_standard   = c() # Standard DDG model
  # Initialize overall net diversification rate vector
  div_rate_log        = c() # Log DDG model
  div_rate_standard   = c() # Standard DDG model
  #
  for (i in 1:length(n_i)){
    spec_rate_log[i]      = rho_w*(log(n_i[i])+1)^{w_d}
    spec_rate_standard[i] = rho_w*(n_i[i]^{w_d})
    #
    ext_rate_log[i]       = rho_e*(log(n_i[i])+1)^{e_d}
    ext_rate_standard[i]  = rho_e*(n_i[i]^{e_d})
    #
    div_rate_log[i]       = spec_rate_log[i] - ext_rate_log[i]
    div_rate_standard[i]  = spec_rate_standard[i] - ext_rate_standard[i]
  }
  # browser()
  #
  # spec_rate_log = log(spec_rate_log)
  # ext_rate_log  = log(ext_rate_log)
  # div_rate_log  = log(div_rate_log)
  # #
  # spec_rate_standard = log(spec_rate_standard)
  # ext_rate_standard  = log(ext_rate_standard)
  # div_rate_standard  = log(div_rate_standard)
  ##############
  # DATA PROCESSING   #
  ##############
  # browser()
  dat_log          <- data.frame(model_type = rep("log",3*length(n_i)),num_species = rep(n_i,3), rates = c(spec_rate_log,ext_rate_log,div_rate_log),
                                 group = rep(c('Speciation rate','Extinction rate',"Net diversification"),each=length(n_i)))
  #
  dat_standard      <- data.frame(model_type = rep("standard",3*length(n_i)),num_species = rep(n_i,3), rates = c(spec_rate_standard,ext_rate_standard,div_rate_standard),
                                  group = rep(c('Speciation rate','Extinction rate',"Net diversification"),each=length(n_i)))
  #
  dat_comb = rbind(dat_log,dat_standard)
  ##############
  # PLOTTINGS  #
  ##############
  #
  title_list <- c("logarithmic DDG","standard DDG")
  # make plot
  p_log = ggplot(dat_log, aes(x = num_species, y = rates, color = group, group = group)) +
          geom_line() +
          labs(title = title_list[1], x = "number of species", y = "rate values")
  #
  p_standard = ggplot(dat_standard, aes(x = num_species, y = rates, color = group, group = group)) +
               geom_line() +
               labs(title = title_list[2], x = "number of species", y = "rate values")
  #
  p_comb = grid.arrange(p_log, p_standard, ncol = 2)  #
  #
  print(p_comb)
  return(dat_comb)
}