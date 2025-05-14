# this script gives a visual representation of distribution of origin age (and root age) and extant tree sizes from our
# training data, and compare it to empirical datasets (viburnum and anolis)

# [speciation_status,extinction_status,dispersal_status]
# 1 = [0,0,0] #geosse
# 2 = [1,0,0] #ddwithin
# 3 = [0,1,0] #ddext
# 4 = [0,0,1] #dddisp
# 5 = [1,1,0] #ddwithin_ext
# 6 = [0,1,1] #ddext_disp
# 7 = [1,0,1] #ddwithin_disp
# 8 = [1,1,1] #ddall
base_path      = Sys.getenv("HOME") # get to user's base directory

in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7_mixture","combined_training_data_onehot")
out_fp         = file.path(base_path,"Downloads/")
setwd(in_fp)

dat        = data.frame(root_age = NA, num_taxa = NA, model_type = NA) #simulated data with root age (x-axis) and num extant taxa (y-axis)

num_sim = 400000

viburnum_dat <- data.frame(root_age = 10.715, num_taxa = 38)  # viburnum's root age and taxa count
anolis_dat <- data.frame(root_age = 64.89, num_taxa = 158) # anolis's root age and taxa count
#add modified viburnum here too e.g., root_age * 5

# fill dat dataframe
for (i in 1:num_sim){
  tre   = read.tree(paste0("sim.",i,".extant.tre")) # load extant tree i
  label =  read.csv(paste0("sim.",i,".labels.csv")) # load char matrix for extant tree i
  #
  get_age = max(branching.times(tre)) # get root age for extant tree i (not origin age)
  get_taxa = length(tre$tip.label) # get number of extant taxa for tree i
  if (label$speciation_status == 0 & label$extinction_status == 0 & label$dispersal_status == 0){
    model_type = "geosse"
  } else if (label$speciation_status == 1 & label$extinction_status == 0 & label$dispersal_status == 0){
    model_type = "ddwithin"
  } else if (label$speciation_status == 0 & label$extinction_status == 1 & label$dispersal_status == 0){
    model_type = "ddext"
  } else if (label$speciation_status == 0 & label$extinction_status == 0 & label$dispersal_status == 1){
    model_type = "dddisp"
  } else if (label$speciation_status == 1 & label$extinction_status == 1 & label$dispersal_status == 0){
    model_type = "ddwithin_ext"
  } else if (label$speciation_status == 0 & label$extinction_status == 1 & label$dispersal_status == 1){
    model_type = "ddext_disp"
  } else if (label$speciation_status == 1 & label$extinction_status == 0 & label$dispersal_status == 1){
    model_type = "ddwithin_disp"
  } else if (label$speciation_status == 1 & label$extinction_status == 1 & label$dispersal_status == 1){
    model_type = "ddall"
  } 
  dat[i,] = c(get_age,get_taxa,model_type) #append to data frame
}

dat$root_age = as.numeric(dat$root_age) # change to numeric
dat$num_taxa = as.numeric(dat$num_taxa) # change to numeric

# fill dat_origin dataframe
dat_origin = dat
names(dat_origin) = c("origin_age","num_taxa","model_type")
for (i in 1:num_sim){
  tre = read.tree(paste0("sim.",i,".extant.tre")) # load extant tree i
  origin_age = dat_origin$origin_age[i] + tre$root.edge #replace root age for tree i with its origin age 
  dat_origin$origin_age[i] = origin_age
}

dat_geosse = dat[dat$model_type == "geosse", ]
dat_within = dat[dat$model_type == "ddwithin", ]
dat_ext = dat[dat$model_type == "ddext", ]
dat_disp = dat[dat$model_type == "dddisp", ]
dat_within_ext = dat[dat$model_type == "ddwithin_ext", ]
dat_ext_disp = dat[dat$model_type == "ddext_disp", ]
dat_within_disp = dat[dat$model_type == "ddwithin_disp", ]
dat_all = dat[dat$model_type == "ddall", ]

p_geosse = ggplot(dat_geosse, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 

p_within = ggplot(dat_within, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 

p_ext = ggplot(dat_ext, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 

p_disp = ggplot(dat_disp, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 

p_within_ext = ggplot(dat_within_ext, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 

p_ext_disp = ggplot(dat_ext_disp, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 

p_within_disp = ggplot(dat_within_disp, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 

p_all = ggplot(dat_all, aes(x = root_age, y = num_taxa)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "red", size = 3, shape = 16) +
  theme_minimal() 


p_combined = ggplot(dat, aes(x = root_age, y = num_taxa,color=model_type)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "black", size = 3, shape = 16) +
  theme_minimal() 

p_combined_origin = ggplot(dat_origin, aes(x = origin_age, y = num_taxa,color=model_type)) +
  geom_point() +
  geom_point(data = viburnum_dat, aes(x = root_age, y = num_taxa), color = "blue", size = 3, shape = 16) +
  geom_point(data = anolis_dat, aes(x = root_age, y = num_taxa), color = "black", size = 3, shape = 16) +
  theme_minimal() 

ggsave(paste0(out_fp,"geosse_scatter.png"), p_geosse, height = 10, width = 15)
ggsave(paste0(out_fp,"within_scatter.png"), p_within, height = 10, width = 15)
ggsave(paste0(out_fp,"ext_scatter.png"), p_ext, height = 10, width = 15)
ggsave(paste0(out_fp,"disp_scatter.png"), p_disp, height = 10, width = 15)
ggsave(paste0(out_fp,"within_ext_scatter.png"), p_within_ext, height = 10, width = 15)
ggsave(paste0(out_fp,"ext_disp_scatter.png"), p_ext_disp, height = 10, width = 15)
ggsave(paste0(out_fp,"within_disp_scatter.png"), p_within_disp, height = 10, width = 15)
ggsave(paste0(out_fp,"ddall.png"), p_all, height = 10, width = 15)

ggsave(paste0(out_fp,"combined.png"), p_combined, height = 10, width = 15)
ggsave(paste0(out_fp,"combined_origin.png"), p_combined_origin, height = 10, width = 15)