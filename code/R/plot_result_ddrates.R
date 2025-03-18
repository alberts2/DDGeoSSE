# Show trajectory under the case of no other factors except density-dependent, and each region share
# equal local equiilibrium capacity (independent of which region)
# Note: 
# The x-axis in the plot shows the number of species change in current (away) region j
# For sum of dispersal into j: we assume all the other regions also have the same num of species as num_spec_home
#
# num_reg : number of regions in the system
# type: "simulate" if it is true generating values "estimate" if it is estimated from phyddle analyses
# folder_name: the main folder where all the other directories created by phyddle analysis are located
# file_name: name of file being read (only for "individual" summary_type)
# min_species_curr : lower bound on number of species occupying a current region j
# max_species_curr: upper bound on number of species occupying a current region j
# num_spec_home: number of species in home region i
# summary_type: "batch" (take average across all trees and plot) or "individual" (over an individual tree)

#e.g.
# plot_rates_flow(5,"estimate","ddall_2","sim.empirical_est.labels_num.csv",100,1000,1,summary_type = "individual")

plot_rates_flow <- function(num_reg,type,folder_name,file_name = NULL,min_species_curr,max_species_curr,num_spec_home,summary_type){
  ##############
  # LOAD R LIBRARIES #
  ##############
  library(gridExtra)
  library(ggplot2)
  #
  ##############
  # FILE SYSTEM #
  ##############
  # Check if file_name and folder_name are in appropriate type (chr)
  if (is.character(folder_name) == FALSE){
    warning("folder_name has to be a chr type")
  } else {
    
  }
  #
  if (is.character(file_name) == FALSE){
    warning("folder_name has to be a chr type")
  } else {
    
  }
  # browser()
  # NOTE: when making this as a package, the design will just be home/folder_name
  base_path      = Sys.getenv("HOME") # get to user's base directory
  fp             = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_3",folder_name)
  if (file.exists(fp) == TRUE){
    in_fp          = paste0(fp,'/',type,'/')
  } else {
    warning(paste0("path ",fp, " does not exist"))
  }
  # browser()
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
          #### Get rate scalar for outgoing dispersal d_i
          if ("exp_scalar_d_i_value" %in% names(param_dat) == TRUE){
            d_i = log(param_dat[1,"exp_scalar_d_i_value"]) # rate scalar for outgoing dispersal d_i
          } else {
            d_i = 0 # must be the case where no density effect on d_i process 
          }
          #### Get rate scalar for incoming dispersal d_j
          if ("exp_scalar_d_j_value" %in% names(param_dat) == TRUE){
            d_j = log(param_dat[1,"exp_scalar_d_j_value"]) # rate scalar for outgoing dispersal d_j
          } else {
            d_j = 0 # must be the case where no density effect on d_j process 
          }
          #### Get rate scalar for between-region speciation b
          if ("exp_scalar_b_value" %in% names(param_dat) == TRUE){
            b_d = log(param_dat[1,"exp_scalar_b_value"]) # rate scalar for between-region b
          } else {
            b_d = 0 # must be the case where no density effect on b process
          }
          #### Get base rate rho_w
          if ("log_base_w_value" %in% names(param_dat) == TRUE){
            rho_w = exp(param_dat[1,"log_base_w_value"]) # base rate for within speciation
          } else {
            warning("Invalid data, does not contain within-region speciation process")
          }
          #### Get base rate rho_e
          if ("log_base_e_value" %in% names(param_dat) == TRUE){
            rho_e = exp(param_dat[1,"log_base_e_value"]) # base rate for extinction
          } else {
            warning("Invalid data, does not contain extinction process")
          }
          #### Get base rate rho_d
          if ("log_base_d_value" %in% names(param_dat) == TRUE){
            rho_d = exp(param_dat[1,"log_base_d_value"]) # base rate for dispersal
          } else {
            warning("Invalid data, does not contain dispersal process")
          }
          #### Get base rate rho_b
          if ("log_base_b_value" %in% names(param_dat) == TRUE){
            rho_b = exp(param_dat[1,"log_base_b_value"]) # base rate for between-region speciation
          } else {
            warning("Invalid data, does not contain between-region speciation process")
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
          #### Get rate scalar for d_i
          if ("exp_scalar_d_i" %in% names(param_dat) == TRUE){
            d_i = log(param_dat[1,"exp_scalar_d_i"]) # rate scalar for outgoing dispersal d_i
          } else {
            warning("invalid data, does not contain rate scalar for outgoing dispersal")
          }
          #### Get rate scalar for d_j
          if ("exp_scalar_d_j" %in% names(param_dat) == TRUE){
            d_j = log(param_dat[1,"exp_scalar_d_j"]) # rate scalar for incoming dispersal d_j
          } else {
            warning("invalid data, does not contain rate scalar for incoming dispersal")
          }
          #### Get rate scalar for b
          if ("exp_scalar_b" %in% names(param_dat) == TRUE){
            b_d = log(param_dat[1,"exp_scalar_b"]) # rate scalar for between-region b
          } else {
            warning("invalid data, does not contain rate scalar for between-region speciation")
          }
          #### Get base rate rho_w
          if ("log_base_w" %in% names(param_dat) == TRUE){
            rho_w = exp(param_dat[1,"log_base_w"]) # base rate for within speciation
          } else {
            warning("Invalid data, does not contain within-region speciation process")
          }
          #### Get base rate rho_e
          if ("log_base_e" %in% names(param_dat) == TRUE){
            rho_e = exp(param_dat[1,"log_base_e"]) # base rate for extinction
          } else {
            warning("Invalid data, does not contain extinction process")
          }
          #### Get base rate rho_d
          if ("log_base_d" %in% names(param_dat) == TRUE){
            rho_d = exp(param_dat[1,"log_base_d"]) # base rate for dispersal
          } else {
            warning("Invalid data, does not contain dispersal process")
          }
          #### Get base rate rho_b
          if ("log_base_b" %in% names(param_dat) == TRUE){
            rho_b = exp(param_dat[1,"log_base_b"]) # base rate for between-region speciation
          } else {
            warning("Invalid data, does not contain between-region speciation process")
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
          d_i_vec     = c() # vector of d_i across all simulations
          d_j_vec     = c() # vector of d_j across all simulations 
          b_d_vec     = c() # vector of b across all simulations
          rho_w_vec   = c() # vector of rho_w across all simulations
          rho_e_vec   = c() # vector of rho_e across all simulations 
          rho_d_vec   = c() # vector of rho_d across all simulations 
          rho_b_vec   = c() # vector of rho_b across all simulations 
          # loop over these files
          for (i in length(indices_vec)){
            # read param data for sim i
            param_dat = read.csv(paste0("sim.",indices_vec[i],".labels.csv"))
            #
            # Get rate scalar for w
            if ("exp_scalar_w" %in% names(param_dat) == TRUE){
              w_d = log(param_dat[1,"exp_scalar_w"]) # rate scalar for within speciation
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain rate scalar for within speciation"))
            }
            w_d_vec = append(w_d_vec,w_d) # append w_d from sim_i to w_d_vec
            # Get rate scalar for e
            if ("exp_scalar_e" %in% names(param_dat) == TRUE){
              e_d = log(param_dat[1,"exp_scalar_e"]) # rate scalar for extinction
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain rate scalar for extinction"))
            }
            e_d_vec = append(e_d_vec,e_d) # append e_d from sim_i to e_d_vec
            #### Get rate scalar for d_i
            if ("exp_scalar_d_i" %in% names(param_dat) == TRUE){
              d_i = log(param_dat[1,"exp_scalar_d_i"]) # rate scalar for outgoing dispersal d_i
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain rate scalar for outgoing dispersal"))
            }
            d_i_vec = append(d_i_vec,d_i) # append d_i from sim_i to d_i_vec
            #### Get rate scalar for d_j
            if ("exp_scalar_d_j" %in% names(param_dat) == TRUE){
              d_j = log(param_dat[1,"exp_scalar_d_j"]) # rate scalar for incoming dispersal d_j
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain rate scalar for incoming dispersal"))
            }
            d_j_vec = append(d_j_vec,d_j) # append d_j from sim_i to d_j_vec
            #### Get rate scalar for b
            if ("exp_scalar_b" %in% names(param_dat) == TRUE){
              b_d = log(param_dat[1,"exp_scalar_b"]) # rate scalar for between-region b
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain rate scalar for between-region speciation"))
            }
            b_d_vec = append(b_d_vec,b_d) # append b_d from sim_i to b_d_vec
            # Get base rate rho_w
            if ("log_base_w" %in% names(param_dat) == TRUE){
              rho_w = exp(param_dat[1,"log_base_w"]) # base rate for within speciation
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain within-region speciation process"))
            }
            rho_w_vec = append(rho_w_vec,rho_w)
            # Get base rate rho_e
            if ("log_base_e" %in% names(param_dat) == FALSE){
              rho_e = exp(param_dat[1,"log_base_e"]) # base rate for extinction
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain extinction process"))
            }
            rho_e_vec = append(rho_e_vec,rho_e)
            #### Get base rate rho_d
            if ("log_base_d" %in% names(param_dat) == TRUE){
              rho_d = exp(param_dat[1,"log_base_d"]) # base rate for dispersal
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain dispersal process"))
            }
            rho_d_vec = append(rho_d_vec,rho_d)
            #### Get base rate rho_b
            if ("log_base_b" %in% names(param_dat) == TRUE){
              rho_b = exp(param_dat[1,"log_base_b"]) # base rate for between-region speciation
            } else {
              warning(paste0("file sim.",indices_vec[i],".labels.csv does not contain between-region speciation process"))
            }
            rho_b_vec = append(rho_b_vec,rho_b)
          }
          # get the average
          w_d   = mean(w_d_vec)
          e_d   = mean(e_d_vec)
          d_i   = mean(d_i_vec)
          d_j   = mean(d_j_vec)
          #
          rho_w = mean(rho_w_vec)
          rho_e = mean(rho_e_vec)
          rho_d = mean(rho_d_vec)
          rho_b = mean(rho_b_vec)
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
  if (min_species_curr >= max_species_curr){
    warning("Max species count must be larger than min species count")
  } else {
 
  }
  #
  if (num_spec_home <= 0){
    warning("species count in home region cannot be 0 or negative")
  } else {

  }
  #
  if (min_species_curr <= 0 || max_species_curr <=0){
    warning("Minimum number of species cannot be 0 or negative")
  } else {
    n_j_curr  = seq(min_species_curr,max_species_curr)
  }
  # Initialize overall speciation rate vector in curr region j
  spec_rate_log       = c() # Log DDG model
  spec_rate_standard  = c() # Standard DDG model
  # Initialize overall extinction rate vector in curr region j
  ext_rate_log        = c() # Log DDG model
  ext_rate_standard   = c() # Standard DDG model
  # Initialize overall dispersal rate vector into curr region j
  disp_rate_log        = c() # Log DDG model
  disp_rate_standard   = c() # Standard DDG model
  # Initialize overall net diversification rate vector in j (influx - outflux from local equilibrium equation)
  div_rate_log        = c() # Log DDG model
  div_rate_standard   = c() # Standard DDG model
  #
  for (j in 1:length(n_j_curr)){
    # compute overall speciation rate in curr region j for both models
    spec_rate_log[j]      = rho_w*(log(n_j_curr[j])+1)^{w_d}
    spec_rate_standard[j] = rho_w*(n_j_curr[j]^{w_d})
    # compute overall extinction rate in curr region j for both models
    ext_rate_log[j]       = rho_e*(log(n_j_curr[j])+1)^{e_d}
    ext_rate_standard[j]  = rho_e*(n_j_curr[j]^{e_d})
    #
    disp_rate_log[j]      = rho_d*(log(num_spec_home)+1)^{d_i}*(log(n_j_curr[j]+1)+1)^{d_j}
    disp_rate_standard[j] = rho_d*(num_spec_home)^{d_i}*(n_j_curr[j]+1)^{d_j}
    #
    div_rate_log[j]       = spec_rate_log[j] - ext_rate_log[j] + (num_reg-1)*disp_rate_log[j]
    div_rate_standard[j]  = spec_rate_standard[j] - ext_rate_standard[j] + (num_reg-1)*disp_rate_standard[j]
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
  dat_log          <- data.frame(model_type = rep("log",3*length(n_j_curr)),num_species = rep(n_j_curr,3), rates = c(log(spec_rate_log),log(ext_rate_log),log(disp_rate_log)),
                                 group = rep(c('Speciation rate','Extinction rate','Dispersal rate'),each=length(n_j_curr)))
  
  dat_netflow_log          <- data.frame(model_type = rep("log",1*length(n_j_curr)),num_species = rep(n_j_curr,1), rates = c(div_rate_log),
                                 group = rep(c('Net flow'),each=length(n_j_curr)))
  #
  dat_standard      <- data.frame(model_type = rep("standard",3*length(n_j_curr)),num_species = rep(n_j_curr,3), rates = c(log(spec_rate_standard),log(ext_rate_standard),log(disp_rate_standard)),
                                  group = rep(c('Speciation rate','Extinction rate','Dispersal rate'),each=length(n_j_curr)))
  
  dat_netflow_standard      <- data.frame(model_type = rep("standard",1*length(n_j_curr)),num_species = rep(n_j_curr,1), rates = c(div_rate_standard),
                                  group = rep(c('Net flow'),each=length(n_j_curr)))
  #
  dat_comb = rbind(dat_log,dat_standard)
  ##############
  # PLOTTINGS  #
  ##############
  #
  title_list <- c("logarithmic DDG","standard DDG","net flow logarithmic DDG","net flow standard DDG")
  # make plot
  p_log              = ggplot(dat_log, aes(x = num_species, y = rates, color = group, group = group)) +
                       geom_line() +
                       labs(title = title_list[1], x = "number of species", y = "log rate values", color = "rate types")
  #
  p_standard         = ggplot(dat_standard, aes(x = num_species, y = rates, color = group, group = group)) +
                       geom_line() +
                       labs(title = title_list[2], x = "number of species", y = "log rate values", color = "rate types")
  #
  p_netflow_log      = ggplot(dat_netflow_log, aes(x = num_species, y = rates, color = group, group = group)) +
                       geom_line() +
                       labs(title = title_list[3], x = "number of species", y = "rate values", color = "rate type")
  #
  p_netflow_standard = ggplot(dat_netflow_standard, aes(x = num_species, y = rates, color = group, group = group)) +
                       geom_line() + 
                       labs(title = title_list[4], x = "number of species", y = "rate values", color = "rate type")
  #
  p_comb = grid.arrange(p_log, p_standard,p_netflow_log,p_netflow_standard,nrow = 2, ncol = 2)  #
  #
  print(p_comb)
  return(dat_comb)
}