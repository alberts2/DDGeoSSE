library(ape)
library(ggplot2)

base_path      = Sys.getenv("HOME") # get to user's base directory
# in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_7","ddall_7","simulate")
in_fp          = file.path(base_path,"Documents","Code Testing","density_geosse","data","deep_learning","server","DD-GeoSSE","5_reg","CURRENT_SET","model_selection_2","ddext_2","simulate")

setwd(in_fp)

num_sim = 50000

#### extinction
low_rhoe_low_ed   = 0
low_rhoe_high_ed  = 0
high_rhoe_low_ed  = 0
high_rhoe_high_ed = 0
#
mod_rhoe_low_ed   = 0
mod_rhoe_high_ed  = 0 
low_rhoe_mod_ed   = 0
high_rhoe_mod_ed  = 0
mod_rhoe_mod_ed   = 0
# treesize of extant trees 
low_rhoe_low_ed_treesize   = c()
low_rhoe_high_ed_treesize  = c()
high_rhoe_low_ed_treesize  = c()
high_rhoe_high_ed_treesize = c()
#
mod_rhoe_low_ed_treesize   = c()
mod_rhoe_high_ed_treesize  = c()
low_rhoe_mod_ed_treesize   = c()
high_rhoe_mod_ed_treesize  = c()
mod_rhoe_mod_ed_treesize   = c()
# treesize of complete trees 
low_rhoe_low_ed_treesize_complete   = c()
low_rhoe_high_ed_treesize_complete  = c()
high_rhoe_low_ed_treesize_complete  = c()
high_rhoe_high_ed_treesize_complete = c()
#
mod_rhoe_low_ed_treesize_complete   = c()
mod_rhoe_high_ed_treesize_complete  = c()
low_rhoe_mod_ed_treesize_complete   = c()
high_rhoe_mod_ed_treesize_complete  = c()
mod_rhoe_mod_ed_treesize_complete   = c()
#### 
#
# high_ed_vec = 0

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
  dat_complete_tree = read.tree(paste0("sim.",i,".tre"))
  dat_extant_tree   = read.tree(paste0("sim.",i,".extant.tre"))
  #
  if (log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)<low_rhoe){ 
    low_rhoe_low_ed = low_rhoe_low_ed + 1
    low_rhoe_low_ed_treesize = append(low_rhoe_low_ed_treesize,length(dat_extant_tree$tip.label))
    low_rhoe_low_ed_treesize_complete = append(low_rhoe_low_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)>high_rhoe){
    high_rhoe_low_ed = high_rhoe_low_ed + 1
    high_rhoe_low_ed_treesize = append(high_rhoe_low_ed_treesize,length(dat_extant_tree$tip.label))
    high_rhoe_low_ed_treesize_complete = append(high_rhoe_low_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<low_rhoe){
    low_rhoe_high_ed = low_rhoe_high_ed + 1
    low_rhoe_high_ed_treesize = append(low_rhoe_high_ed_treesize,length(dat_extant_tree$tip.label))
    low_rhoe_high_ed_treesize_complete = append(low_rhoe_high_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)>high_rhoe){
    high_rhoe_high_ed = high_rhoe_high_ed + 1
    high_rhoe_high_ed_treesize = append(high_rhoe_high_ed_treesize,length(dat_extant_tree$tip.label))
    high_rhoe_high_ed_treesize_complete = append(high_rhoe_high_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe){
    mod_rhoe_low_ed = mod_rhoe_low_ed + 1
    mod_rhoe_low_ed_treesize = append(mod_rhoe_low_ed_treesize,length(dat_extant_tree$tip.label))
    mod_rhoe_low_ed_treesize_complete = append(mod_rhoe_low_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) > high_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe){
    mod_rhoe_high_ed = mod_rhoe_high_ed + 1
    mod_rhoe_high_ed_treesize = append(mod_rhoe_high_ed_treesize,length(dat_extant_tree$tip.label))
    mod_rhoe_high_ed_treesize_complete = append(mod_rhoe_high_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)<low_rhoe){
    low_rhoe_mod_ed = low_rhoe_mod_ed + 1
    low_rhoe_mod_ed_treesize = append(low_rhoe_mod_ed_treesize,length(dat_extant_tree$tip.label))
    low_rhoe_mod_ed_treesize_complete = append(low_rhoe_mod_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)>high_rhoe){
    high_rhoe_mod_ed = high_rhoe_mod_ed + 1
    high_rhoe_mod_ed_treesize = append(high_rhoe_mod_ed_treesize,length(dat_extant_tree$tip.label))
    high_rhoe_mod_ed_treesize_complete = append(high_rhoe_mod_ed_treesize_complete,length(dat_complete_tree$tip.label))
  } else if (log(dat_labels$exp_scalar_e) < high_ed && log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)<high_rhoe && exp(dat_labels$log_base_e)>low_rhoe){
    mod_rhoe_mod_ed = mod_rhoe_mod_ed+1
    mod_rhoe_mod_ed_treesize = append(mod_rhoe_mod_ed_treesize,length(dat_extant_tree$tip.label))
    mod_rhoe_mod_ed_treesize_complete = append(mod_rhoe_mod_ed_treesize_complete,length(dat_complete_tree$tip.label))
  }
}

# for (i in 1:num_sim){
#   dat_labels        = read.csv(paste0("sim.",i,".labels.csv")) # read labels file
#   #
#   if (log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)<low_rhoe){ 
#     low_rhoe_low_ed = low_rhoe_low_ed + 1
#   } else if (log(dat_labels$exp_scalar_e) < low_ed && exp(dat_labels$log_base_e)>low_rhoe){
#     high_rhoe_low_ed = high_rhoe_low_ed + 1
#   } else if (log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)<low_rhoe){
#     low_rhoe_high_ed = low_rhoe_high_ed + 1
#   } else if (log(dat_labels$exp_scalar_e) > low_ed && exp(dat_labels$log_base_e)>low_rhoe){
#     high_rhoe_high_ed = high_rhoe_high_ed + 1
#   }
# }

# for (i in 1:num_sim){
#   dat_labels        = read.csv(paste0("sim.",i,".labels.csv")) # read labels file
#   if (log(dat_labels$exp_scalar_e) > low_ed){
#     high_ed_vec = high_ed_vec+1
#   }
# }

# Sample data
data_ext <- data.frame(
  category = c("Low_pe_Low_ed", "Low_pe_High_ed", "High_pe_Low_ed", "High_pe_High_ed",
               "Mod_pe_Low_ed", "Mod_pe_High_ed","Low_pe_Mod_ed","High_pe_Mod_ed","Mod_pe_Mod_ed"),
  value = c(low_rhoe_low_ed, low_rhoe_high_ed, high_rhoe_low_ed, high_rhoe_high_ed,
            mod_rhoe_low_ed,mod_rhoe_high_ed,low_rhoe_mod_ed,high_rhoe_mod_ed,mod_rhoe_mod_ed),
  taxa  = c(mean(low_rhoe_low_ed_treesize),mean(low_rhoe_high_ed_treesize),mean(high_rhoe_low_ed_treesize),mean(high_rhoe_high_ed_treesize),
            mean(mod_rhoe_low_ed_treesize),mean(mod_rhoe_high_ed_treesize),mean(low_rhoe_mod_ed_treesize),mean(high_rhoe_mod_ed_treesize),
            mean(mod_rhoe_mod_ed_treesize)),
  taxa_loss = c(mean((low_rhoe_low_ed_treesize_complete-low_rhoe_low_ed_treesize)/low_rhoe_low_ed_treesize_complete)*100,
                mean((low_rhoe_high_ed_treesize_complete-low_rhoe_high_ed_treesize)/low_rhoe_high_ed_treesize_complete)*100,
                mean((high_rhoe_low_ed_treesize_complete-high_rhoe_low_ed_treesize)/high_rhoe_low_ed_treesize_complete)*100,
                mean((high_rhoe_high_ed_treesize_complete-high_rhoe_high_ed_treesize)/high_rhoe_high_ed_treesize_complete)*100,
                mean((mod_rhoe_low_ed_treesize_complete-mod_rhoe_low_ed_treesize)/mod_rhoe_low_ed_treesize_complete)*100,
                mean((mod_rhoe_high_ed_treesize_complete-mod_rhoe_high_ed_treesize)/mod_rhoe_high_ed_treesize_complete)*100,
                mean((low_rhoe_mod_ed_treesize_complete-low_rhoe_mod_ed_treesize)/low_rhoe_mod_ed_treesize_complete)*100,
                mean((high_rhoe_mod_ed_treesize_complete-high_rhoe_mod_ed_treesize)/high_rhoe_mod_ed_treesize_complete)*100,
                mean((mod_rhoe_mod_ed_treesize_complete-mod_rhoe_mod_ed_treesize)/mod_rhoe_mod_ed_treesize_complete)*100)
)

data_ext$percentage <- round(100 * data_ext$value / sum(data_ext$value), 1) # Calculate percentage
data_ext$label <- paste0(data_ext$category, " (", data_ext$percentage, "%)"," (",round(data_ext$taxa,2),")"," (",round(data_ext$taxa_loss,2),"%)") # Category and percentage in one line

# Compute position for labels
data_ext$ypos <- cumsum(data_ext$value) - data_ext$value / 2


# Create a pie chart using ggplot2
ggplot(data_ext, aes(x = "Extinction Rate", y = value, fill = category)) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) + # Labels inside bars
  theme_minimal() + 
  labs(title = "Extinction Rate Variation", x = NULL, y = "Value") +  
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))