# base_path      = Sys.getenv("HOME") # get to user's base directory
# 
# in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7_mixture","Anoles","cat_estimate_7_anoles","simulate")
# out_fp         = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7_mixture","Anoles","cat_estimate_7_anoles","simulate/")
# 
# setwd(in_fp)
# 
# for (i in 100001:150000){
#   dat_labels        = read.csv(paste0("sim.",i,".labels.csv")) # read labels file
#   #
#   dat_labels$model_type = 1
#   #
#   write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
# }

# This script relabels model type from simulated dataset from integer to one-hot encoding.
# For example
# [speciation_status,extinction_status,dispersal_status]
# 1 = [0,0,0] #geosse
# 2 = [1,0,0] #ddwithin
# 3 = [0,1,0] #ddext
# 4 = [0,0,1] #dddisp
# 5 = [1,1,0] #ddwithin_ext
# 6 = [0,1,1] #ddext_disp
# 7 = [1,0,1] #ddwithin_disp
# 8 = [1,1,1] #ddall
library(ape)

base_path      = Sys.getenv("HOME") # get to user's base directory

in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7_mixture","Anoles","cat_estimate_7_anolis","N1","simulate")
out_fp         = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7_mixture","Anoles","cat_estimate_7_anolis","N3_one_hot","simulate/")
setwd(in_fp)

for (i in 1:400000){
  dat_labels        = read.csv(paste0("sim.",i,".labels.csv")) # read labels file
  dat_char          = read.csv(paste0("sim.",i,".dat.csv")) # read character file
  dat_phy           = read.tree(paste0("sim.",i,".tre")) # read the phylogeny
  #
  if (dat_labels$model_type == 1){
    dat_labels$model_type = NULL
    #
    dat_labels$speciation_status = 0
    dat_labels$extinction_status = 0
    dat_labels$dispersal_status = 0
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),digits = 50) # need to adjust the digits otherwise will get incorrect branch lengths for an extinct taxa really close to the present
  } else if (dat_labels$model_type == 2) {
    dat_labels$model_type = NULL
    #
    dat_labels$speciation_status = 1
    dat_labels$extinction_status = 0
    dat_labels$dispersal_status = 0
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),digits = 50)
  } else if (dat_labels$model_type == 3) {
    dat_labels$model_type = NULL
    #
    dat_labels$speciation_status = 0
    dat_labels$extinction_status = 1
    dat_labels$dispersal_status = 0
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),digits = 50)
  } else if (dat_labels$model_type == 4) {
    dat_labels$model_type = NULL
    #
    dat_labels$speciation_status = 0
    dat_labels$extinction_status = 0
    dat_labels$dispersal_status = 1
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),digits = 50)
  } else if (dat_labels$model_type == 5) {
    dat_labels$model_type = NULL
    #
    dat_labels$speciation_status = 1
    dat_labels$extinction_status = 1
    dat_labels$dispersal_status = 0
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),digits = 50)
  } else if (dat_labels$model_type == 6) {
    dat_labels$model_type = NULL
    #
    dat_labels$speciation_status = 0
    dat_labels$extinction_status = 1
    dat_labels$dispersal_status = 1
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),digits = 50)
  } else if (dat_labels$model_type == 7) {
    dat_labels$model_type = NULL
    #
    dat_labels$speciation_status = 1
    dat_labels$extinction_status = 0
    dat_labels$dispersal_status = 1
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),digits = 50)
  } else if (dat_labels$model_type == 8) {
    dat_labels$model_type = NULL
    #x
    dat_labels$speciation_status = 1
    dat_labels$extinction_status = 1
    dat_labels$dispersal_status = 1
    write.csv(dat_labels,file = paste0(out_fp,"sim.",i,".labels.csv"),row.names = FALSE )
    write.csv(dat_char,file = paste0(out_fp,"sim.",i,".dat.csv"),row.names = FALSE )
    write.tree(dat_phy,file = paste0(out_fp,"sim.",i,".tre"),,digits = 50)
  }
}
