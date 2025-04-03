library(ape)

base_path      = Sys.getenv("HOME") # get to user's base directory
# in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_3","ddall_3","simulate")
# out_fp         = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_3","filtered_ddall_3","simulate/")

# in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7","ddall_7","simulate")
# out_fp         = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7","filtered_ddall_7","simulate/")

in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_2","ddext_2","simulate")
out_fp         = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_2","filtered_ddext_2","simulate/")

setwd(in_fp)

num_sim = 50000

# compute first and third quartiles from uniform distribution and loguniform distribution
# loguniform is used to draw the base rates 
# uniform is used to draw the density scalars 
# low_ed    = 0.0001+(1/4)*(2+0.0001)
# high_ed   = 0.0001+(3/4)*(2+0.0001)
# low_rhoe  = exp(log(0.0001)+(1/4)*(log(1)-log(0.0001)))
# high_rhoe = exp(log(0.0001)+(3/4)*(log(1)-log(0.0001)))

low_ed    = 0.0001+(1/4)*(2+0.0001)
high_ed   = 0.0001+(3/4)*(2+0.0001)
low_rhoe  = exp(log(0.0001)+(1/4)*(log(5)-log(0.0001)))
high_rhoe = exp(log(0.0001)+(3/4)*(log(5)-log(0.0001)))

for (i in 1:num_sim){
  dat_labels        = read.csv(paste0("sim.",i,".labels.csv")) # read labels file
  dat_char          = read.csv(paste0("sim.",i,".dat.csv")) # read character file
  dat_extant_char   = read.csv(paste0("sim.",i,".extant.dat.csv")) # read extant character file
  dat_phy           = read.tree(paste0("sim.",i,".tre")) # read the phylogeny
  # serial_including_modmod_excluding_lowtaxapruned
  # if ((log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe)||(log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe)){
  #     write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
  #     write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
  #     write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"))
  # }
  # serial_notincluding_modmod_excluding_lowtaxapruned
  if ((log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe)){
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"))
  }
  # serial_notincluding_modmod_including_lowtaxapruned
  # if ((log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)<low_rhoe)||(log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<low_rhoe)||(log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)>high_rhoe)||(log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe)||(log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe)||(log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)<low_rhoe)||(log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)>high_rhoe)){
  #   write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
  #   write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
  #   write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"))
  # }
  # nohigh_rhoe_high_ed
  # if ((log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<low_rhoe)||(log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe)){
  #   write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
  #   write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
  #   write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"))
  # }
}