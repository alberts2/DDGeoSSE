library(ape)

base_path      = Sys.getenv("HOME") # get to user's base directory
in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_3","ddall_3","simulate")
out_fp         = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_3","filtered_ddall_3","simulate/")

setwd(in_fp)

num_sim = 50000
  
for (i in 1:num_sim){
  dat_labels        = read.csv(paste0("sim.",i,".labels.csv")) # read labels file
  dat_char          = read.csv(paste0("sim.",i,".dat.csv")) # read character file
  dat_extant_char   = read.csv(paste0("sim.",i,".extant.dat.csv")) # read extant character file
  dat_phy           = read.tree(paste0("sim.",i,".tre")) # read the phylogeny
  #
  if (log(dat_labels$exp_scalar_e) > 2){ #filtered simulation with strong density effect on extinction
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv") )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv") )
    write.csv(dat_extant_char,file = paste0(out_fp,"sim.",i,".extant.dat.csv") )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"))
  }
}