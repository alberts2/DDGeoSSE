# this script rescalling the root age of viburnum tree to be similar to origin age from our training data
# this results in accurate and consistent prediction for each density feature across two independently trained network
# under the same training data and network architecture

base_path      = Sys.getenv("HOME") # get to user's base directory
in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","empiric_viburnum","viburnum_dat_edit")
out_fp         = file.path(base_path,"Documents","Code Testing","density_geosse","data","empiric_viburnum_rescaled/")

setwd(in_fp)

tree_viburnum = read.tree("oreino5.1.tre")
tree_viburnum$edge.length = tree_viburnum$edge.length*7

write.tree(tree_viburnum,file = paste0(out_fp,"oreino5.1.tre"),digits = 50)