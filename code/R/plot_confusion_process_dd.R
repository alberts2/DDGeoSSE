# This script is similar to plot_confusion_submodel.R
# but instead of creating confusion matrix for submodel
# we create confusion matrix for process-level selection
# the idea is: to measure how much of the confusion comes from a particular feature
# effect parameter.
# so we create 3x3 where true process are dataset simulated under pure w^D, e^D or d^D
# and the predicted process is detecting w^D using network trained to infer w^D, detecting
# e^D using network trained to infer e^D, etc. 

library(ggplot2)
library(reshape2)

###### FILE SETTING ######
fp              = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
dat_fp          = paste0(fp,"data/validation_confusion/")
out_fp          = paste0(fp,"plot/confusion_submodel/")

# network 1
wd_fp           = paste0(dat_fp,"network_1/wd_1/")
ed_fp           = paste0(dat_fp,"network_1/ed_1/")
dd_fp           = paste0(dat_fp,"network_1/dd_1/")

# network 2
wd_fp           = paste0(dat_fp,"network_2/wd_2/")
ed_fp           = paste0(dat_fp,"network_2/ed_2/")
dd_fp           = paste0(dat_fp,"network_2/dd_2/")

# network 3
wd_fp           = paste0(dat_fp,"network_3/wd_3/")
ed_fp           = paste0(dat_fp,"network_3/ed_3/")
dd_fp           = paste0(dat_fp,"network_3/dd_3/")

# network 4
wd_fp           = paste0(dat_fp,"network_4/wd_4/")
ed_fp           = paste0(dat_fp,"network_4/ed_4/")
dd_fp           = paste0(dat_fp,"network_4/dd_4/")

# network 5
wd_fp           = paste0(dat_fp,"network_5/wd_5/")
ed_fp           = paste0(dat_fp,"network_5/ed_5/")
dd_fp           = paste0(dat_fp,"network_5/dd_5/")

# network 6
wd_fp           = paste0(dat_fp,"network_6/wd_6/")
ed_fp           = paste0(dat_fp,"network_6/ed_6/")
dd_fp           = paste0(dat_fp,"network_6/dd_6/")

# network 7
wd_fp           = paste0(dat_fp,"network_7/wd_7/")
ed_fp           = paste0(dat_fp,"network_7/ed_7/")
dd_fp           = paste0(dat_fp,"network_7/dd_7/")

# network 8
wd_fp           = paste0(dat_fp,"network_8/wd_8/")
ed_fp           = paste0(dat_fp,"network_8/ed_8/")
dd_fp           = paste0(dat_fp,"network_8/dd_8/")

# network 1 1000 test data (plot look similar to sample size of 100)
wd_fp           = paste0(dat_fp,"network_1_1000/wd_1_1000/")
ed_fp           = paste0(dat_fp,"network_1_1000/ed_1_1000/")
dd_fp           = paste0(dat_fp,"network_1_1000/dd_1_1000/")

# network 1 1000 test data from mixture logU and U  
wd_fp           = paste0(dat_fp,"network_1_1000_mix/wd_1_mix/")
ed_fp           = paste0(dat_fp,"network_1_1000_mix/ed_1_mix/")
dd_fp           = paste0(dat_fp,"network_1_1000_mix/dd_1_mix/")

# network 1 1000 test data from mixture logU and U with testdata from strong DD effect   
wd_fp           = paste0(dat_fp,"network_1_1000_mix_strong/wd_1_mix_strong/")
ed_fp           = paste0(dat_fp,"network_1_1000_mix_strong/ed_1_mix_strong/")
dd_fp           = paste0(dat_fp,"network_1_1000_mix_strong/dd_1_mix_strong/")

# Using train networks based on wd data only, ed data only and dd data only
# instead of combined dataset 
# 1000 test data for confusion
wd_fp           = paste0(dat_fp,"network_retrain_wd_data_only/strong_testdata_1000/")
ed_fp           = paste0(dat_fp,"network_retrain_ed_data_only/strong_testdata_1000/")
dd_fp           = paste0(dat_fp,"network_retrain_dd_data_only/strong_testdata_1000/")

# Using train networks based on wd data only, ed data only and dd data only
# instead of combined dataset with known parameters 
# 1000 test data for confusion
wd_fp           = paste0(dat_fp,"network_retrain_wd_data_only_known_param/strong_testdata_1000/")
ed_fp           = paste0(dat_fp,"network_retrain_ed_data_only_known_param/strong_testdata_1000/")
dd_fp           = paste0(dat_fp,"network_retrain_dd_data_only_known_param/strong_testdata_1000/")

# network 1 1000 test data from mixture logU and U with testdata from strong DD effect with nonoverlapping base rate
# and dd effect 
wd_fp           = paste0(dat_fp,"network_1_mix_strong_better_test_data/wd_1_mix_better_strong/")
ed_fp           = paste0(dat_fp,"network_1_mix_strong_better_test_data/ed_1_mix_better_strong/")
dd_fp           = paste0(dat_fp,"network_1_mix_strong_better_test_data/dd_1_mix_better_strong/")

# Using train networks based on wd data only, ed data only and dd data only
# instead of combined dataset 
# 1000 test data with nonoverlapping base rate and dd effect
wd_fp           = paste0(dat_fp,"network_retrain_wd_data_only_better_test/strong_testdata_1000/")
ed_fp           = paste0(dat_fp,"network_retrain_ed_data_only_better_test/strong_testdata_1000/")
dd_fp           = paste0(dat_fp,"network_retrain_dd_data_only_better_test/strong_testdata_1000/")

# Network 1 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd 
wd_fp           = paste0(dat_fp,"network_1_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_1_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_1_ori_20000/dd_1/")

# Network 2 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd 
wd_fp           = paste0(dat_fp,"network_2_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_2_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_2_ori_20000/dd_1/")

# Network 3 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd 
wd_fp           = paste0(dat_fp,"network_3_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_3_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_3_ori_20000/dd_1/")

# Network 4 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd 
wd_fp           = paste0(dat_fp,"network_4_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_4_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_4_ori_20000/dd_1/")

# Network 5 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd 
wd_fp           = paste0(dat_fp,"network_5_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_5_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_5_ori_20000/dd_1/")

# Network 6 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd
wd_fp           = paste0(dat_fp,"network_6_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_6_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_6_ori_20000/dd_1/")

# Network 7 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd
wd_fp           = paste0(dat_fp,"network_7_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_7_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_7_ori_20000/dd_1/")

# Network 8 using original test dataset from paper (20000)
# there are 3 separate test datasets here, one for network wd, one for ed, and one for dd
wd_fp           = paste0(dat_fp,"network_8_ori_20000/wd_1/")
ed_fp           = paste0(dat_fp,"network_8_ori_20000/ed_1/")
dd_fp           = paste0(dat_fp,"network_8_ori_20000/dd_1/")

# Network 1 using new test dataset with total samples of 20000
wd_fp           = paste0(dat_fp,"network_1_new_20000/wd_1_mix_2500/")
ed_fp           = paste0(dat_fp,"network_1_new_20000/ed_1_mix_2500/")
dd_fp           = paste0(dat_fp,"network_1_new_20000/dd_1_mix_2500/")

# Network 1 using another test dataset taken from the original dataset with 400000
wd_fp           = paste0(dat_fp,"network_1_ori_another_20000/wd_1_ori_2500/")
ed_fp           = paste0(dat_fp,"network_1_ori_another_20000/ed_1_ori_2500/")
dd_fp           = paste0(dat_fp,"network_1_ori_another_20000/dd_1_ori_2500/")

# Network 1 using test data where I fixed bug on same max_taxa and max_time for all trees in the same batch
wd_fp           = paste0(dat_fp,"network_1_new_nobug_20000/wd_1_new_nobug_2500/")
ed_fp           = paste0(dat_fp,"network_1_new_nobug_20000/ed_1_new_nobug_2500/")
dd_fp           = paste0(dat_fp,"network_1_new_nobug_20000/dd_1_new_nobug_2500/")

# Network 1 using test data where I fixed bug on max_taxa/max_time and also random.seed 
wd_fp           = paste0(dat_fp,"network_1_newer_20000/wd_1_newer_2500/")
ed_fp           = paste0(dat_fp,"network_1_newer_20000/ed_1_newer_2500/")
dd_fp           = paste0(dat_fp,"network_1_newer_20000/dd_1_newer_2500/")

# Network 1 using test data where I fixed bug on max_taxa/max_time and also random.seed 
wd_fp           = paste0(dat_fp,"network_1_newer_2_20000/wd_1_newer_2500/")
ed_fp           = paste0(dat_fp,"network_1_newer_2_20000/ed_1_newer_2500/")
dd_fp           = paste0(dat_fp,"network_1_newer_2_20000/dd_1_newer_2500/")

##### SUBMODEL 0
# get true and estimate w_d for submodel 0
true_wd_model_0 = read.table(paste0(wd_fp,"estimate_0/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_0  = read.table(paste0(wd_fp,"estimate_0/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_0 = cbind(true_wd_model_0,est_wd_model_0[,2:3])

# get true and estimate e_d for submodel 0
true_ed_model_0 = read.table(paste0(ed_fp,"estimate_0/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_0  = read.table(paste0(ed_fp,"estimate_0/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_0 = cbind(true_ed_model_0,est_ed_model_0[,2:3])

# get true and estimate d_d for submodel 0
true_dd_model_0 = read.table(paste0(dd_fp,"estimate_0/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_0  = read.table(paste0(dd_fp,"estimate_0/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_0 = cbind(true_dd_model_0,est_dd_model_0[,2:3])

##### SUBMODEL 1
# get true and estimate w_d for submodel 1 
true_wd_model_1 = read.table(paste0(wd_fp,"estimate_1/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_1  = read.table(paste0(wd_fp,"estimate_1/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_1 = cbind(true_wd_model_1,est_wd_model_1[,2:3])

# get true and estimate e_d for submodel 1 
true_ed_model_1 = read.table(paste0(ed_fp,"estimate_1/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_1  = read.table(paste0(ed_fp,"estimate_1/sim.test_est.labels_cat.csv"),sep = ",",header = T)
# est_ed_model_1  = read.table(paste0(ed_fp,"estimate_1/sim.empirical_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_1 = cbind(true_ed_model_1,est_ed_model_1[,2:3])

# get true and estimate d_d for submodel 1 
true_dd_model_1 = read.table(paste0(dd_fp,"estimate_1/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_1  = read.table(paste0(dd_fp,"estimate_1/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_1 = cbind(true_dd_model_1,est_dd_model_1[,2:3])

##### SUBMODEL 2
# get true and estimate w_d for submodel 2 
true_wd_model_2 = read.table(paste0(wd_fp,"estimate_2/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_2  = read.table(paste0(wd_fp,"estimate_2/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_2 = cbind(true_wd_model_2,est_wd_model_2[,2:3])

# get true and estimate e_d for submodel 2 
true_ed_model_2 = read.table(paste0(ed_fp,"estimate_2/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_2  = read.table(paste0(ed_fp,"estimate_2/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_2 = cbind(true_ed_model_2,est_ed_model_2[,2:3])

# get true and estimate d_d for submodel 2 
true_dd_model_2 = read.table(paste0(dd_fp,"estimate_2/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_2  = read.table(paste0(dd_fp,"estimate_2/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_2 = cbind(true_dd_model_2,est_dd_model_2[,2:3])

##### SUBMODEL 3
# get true and estimate w_d for submodel 3 
true_wd_model_3 = read.table(paste0(wd_fp,"estimate_3/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_3  = read.table(paste0(wd_fp,"estimate_3/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_3 = cbind(true_wd_model_3,est_wd_model_3[,2:3])

# get true and estimate e_d for submodel 3 
true_ed_model_3 = read.table(paste0(ed_fp,"estimate_3/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_3  = read.table(paste0(ed_fp,"estimate_3/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_3 = cbind(true_ed_model_3,est_ed_model_3[,2:3])

# get true and estimate d_d for submodel 3 
true_dd_model_3 = read.table(paste0(dd_fp,"estimate_3/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_3  = read.table(paste0(dd_fp,"estimate_3/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_3 = cbind(true_dd_model_3,est_dd_model_3[,2:3])

##### SUBMODEL 4
# get true and estimate w_d for submodel 4 
true_wd_model_4 = read.table(paste0(wd_fp,"estimate_4/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_4  = read.table(paste0(wd_fp,"estimate_4/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_4 = cbind(true_wd_model_4,est_wd_model_4[,2:3])

# get true and estimate e_d for submodel 4 
true_ed_model_4 = read.table(paste0(ed_fp,"estimate_4/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_4  = read.table(paste0(ed_fp,"estimate_4/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_4 = cbind(true_ed_model_4,est_ed_model_4[,2:3])

# get true and estimate d_d for submodel 4 
true_dd_model_4 = read.table(paste0(dd_fp,"estimate_4/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_4  = read.table(paste0(dd_fp,"estimate_4/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_4 = cbind(true_dd_model_4,est_dd_model_4[,2:3])

##### SUBMODEL 5
# get true and estimate w_d for submodel 5
true_wd_model_5 = read.table(paste0(wd_fp,"estimate_5/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_5  = read.table(paste0(wd_fp,"estimate_5/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_5 = cbind(true_wd_model_5,est_wd_model_5[,2:3])

# get true and estimate e_d for submodel 5
true_ed_model_5 = read.table(paste0(ed_fp,"estimate_5/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_5  = read.table(paste0(ed_fp,"estimate_5/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_5 = cbind(true_ed_model_5,est_ed_model_5[,2:3])

# get true and estimate d_d for submodel 5
true_dd_model_5 = read.table(paste0(dd_fp,"estimate_5/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_5  = read.table(paste0(dd_fp,"estimate_5/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_5 = cbind(true_dd_model_5,est_dd_model_5[,2:3])

##### SUBMODEL 6
# get true and estimate w_d for submodel 6
true_wd_model_6 = read.table(paste0(wd_fp,"estimate_6/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_6  = read.table(paste0(wd_fp,"estimate_6/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_6 = cbind(true_wd_model_6,est_wd_model_6[,2:3])

# get true and estimate e_d for submodel 6
true_ed_model_6 = read.table(paste0(ed_fp,"estimate_6/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_6  = read.table(paste0(ed_fp,"estimate_6/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_6 = cbind(true_ed_model_6,est_ed_model_6[,2:3])

# get true and estimate d_d for submodel 6
true_dd_model_6 = read.table(paste0(dd_fp,"estimate_6/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_6  = read.table(paste0(dd_fp,"estimate_6/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_6 = cbind(true_dd_model_6,est_dd_model_6[,2:3])

##### SUBMODEL 7
# get true and estimate w_d for submodel 7
true_wd_model_7 = read.table(paste0(wd_fp,"estimate_7/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_wd_model_7  = read.table(paste0(wd_fp,"estimate_7/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_wd_model_7 = cbind(true_wd_model_7,est_wd_model_7[,2:3])

# get true and estimate e_d for submodel 7
true_ed_model_7 = read.table(paste0(ed_fp,"estimate_7/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_7  = read.table(paste0(ed_fp,"estimate_7/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_7 = cbind(true_ed_model_7,est_ed_model_7[,2:3])

# get true and estimate d_d for submodel 7
true_dd_model_7 = read.table(paste0(dd_fp,"estimate_7/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_dd_model_7  = read.table(paste0(dd_fp,"estimate_7/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_dd_model_7 = cbind(true_dd_model_7,est_dd_model_7[,2:3])

#### ALL MODELS (TEST DATA FROM ALL MODELS). 
# this is using test data combined from all models 
true_ed_model_all = read.table(paste0(ed_fp,"estimate_allmodels/sim.test_true.labels_cat.csv"),sep = ",",header = T)
est_ed_model_all  = read.table(paste0(ed_fp,"estimate_allmodels/sim.test_est.labels_cat.csv"),sep = ",",header = T)
comb_ed_model_all = cbind(true_ed_model_all,est_ed_model_all[,2:3])



# 2x2 confusion matrix using test data from all models 
true <- comb_ed_model_all$extinction_status

# Predicted labels:
# choose the class with the highest predicted probability
pred <- ifelse(comb_ed_model_all$extinction_status_1 > 0.5, 1, 0)

# 2 x 2 confusion matrix
conf_mat <- table(
  Predicted = pred,
  True = true
)
conf_prop = prop.table(conf_mat, margin = 2)
# Convert to data frame for ggplot
plot_dat <- melt(conf_prop)
colnames(plot_dat) <- c("Predicted", "True", "Probability")



# using mean of count on how many samples with DD turned on or off (>0.5 treshold)
confuse_dat = data.frame(
  true_process = c("none","w_only", "e_only", "d_only","we_only","ed_only","wd_only","all"),
  w_network = c(
    mean(comb_wd_model_0$speciation_status_1>0.5),
    mean(comb_wd_model_2$speciation_status_1>0.5),
    mean(comb_wd_model_1$speciation_status_1>0.5),
    mean(comb_wd_model_3$speciation_status_1>0.5),
    mean(comb_wd_model_4$speciation_status_1>0.5),
    mean(comb_wd_model_5$speciation_status_1>0.5),
    mean(comb_wd_model_6$speciation_status_1>0.5),
    mean(comb_wd_model_7$speciation_status_1>0.5)),
  e_network = c(
    mean(comb_ed_model_0$extinction_status_1>0.5),
    mean(comb_ed_model_2$extinction_status_1>0.5),
    mean(comb_ed_model_1$extinction_status_1>0.5),
    mean(comb_ed_model_3$extinction_status_1>0.5),
    mean(comb_ed_model_4$extinction_status_1>0.5),
    mean(comb_ed_model_5$extinction_status_1>0.5),
    mean(comb_ed_model_6$extinction_status_1>0.5),
    mean(comb_ed_model_7$extinction_status_1>0.5)),
  d_network = c( 
    mean(comb_dd_model_0$dispersal_status_1>0.5),
    mean(comb_dd_model_2$dispersal_status_1>0.5),
    mean(comb_dd_model_1$dispersal_status_1>0.5),
    mean(comb_dd_model_3$dispersal_status_1>0.5),
    mean(comb_dd_model_4$dispersal_status_1>0.5),
    mean(comb_dd_model_5$dispersal_status_1>0.5),
    mean(comb_dd_model_6$dispersal_status_1>0.5),
    mean(comb_dd_model_7$dispersal_status_1>0.5))
)
############
process_long <- data.frame(
  true_process = rep(confuse_dat$true_process, times = 3),
  predicted_process = rep(c("w_network","e_network","d_network"),
                          each = nrow(confuse_dat)),
  prob = c(confuse_dat$w_network,
           confuse_dat$e_network,
           confuse_dat$d_network)
)

process_long$true_process <- factor(
  process_long$true_process,
  levels = c("none","w_only", "e_only", "d_only","we_only","ed_only","wd_only","all")
)

process_long$predicted_process <- factor(
  process_long$predicted_process,
  levels = c("w_network","e_network","d_network")
)

label_fun_process <- function(x) {
  labs <- c(
    none    = "w^D == 0 ~ ',' ~ e^D == 0 ~ ',' ~ d^{list(D, dest)} == 0",
    w_only  = "w^D < 0 ~ ',' ~ e^D == 0 ~ ',' ~ d^{list(D, dest)} == 0",
    e_only  = "w^D == 0 ~ ',' ~ e^D > 0 ~ ',' ~ d^{list(D, dest)} == 0",
    d_only  = "w^D == 0 ~ ',' ~ e^D == 0 ~ ',' ~ d^{list(D, dest)} < 0",
    we_only = "w^D < 0 ~ ',' ~ e^D > 0 ~ ',' ~ d^{list(D, dest)} == 0",
    ed_only = "w^D == 0 ~ ',' ~ e^D > 0 ~ ',' ~ d^{list(D, dest)} < 0",
    wd_only = "w^D < 0 ~ ',' ~ e^D == 0 ~ ',' ~ d^{list(D, dest)} < 0",
    all     = "w^D < 0 ~ ',' ~ e^D > 0 ~ ',' ~ d^{list(D, dest)} < 0",
    w_network = "Detect~w^D",
    e_network = "Detect~e^D",
    d_network = "Detect ~ d^{list(D, dest)}"
  )
  parse(text = labs[x])
}
###########

plot_2by2 = ggplot(plot_dat,
                   aes(x = factor(True),
                       y = factor(Predicted),
                       fill = Probability)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Probability)),
            size = 8) +
  scale_fill_gradient(
    low = "white",
    high = "#8e239b",
    limits = c(0, 1)
  ) +
  labs(
    title = "Test estimates: extinction_status",
    x = "extinction_status truth",
    y = "extinction_status estimate"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

plot_p = ggplot(process_long,
                aes(x = predicted_process,
                    y = true_process,
                    fill = prob)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", prob)), size = 5) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1)) +
  scale_x_discrete(labels = label_fun_process) +
  scale_y_discrete(labels = label_fun_process) +
  theme_minimal() +
  labs(x = "Predicted process",
       y = "True process",
       fill = "Detection probability") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16)
  )


# pdf(file = paste0(out_fp,"confusion_network_1_process_mixture_strong.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_1_process_mixture_strong_better_testdata.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_1_process_mixture.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_1_ori_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_2_ori_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_3_ori_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_4_ori_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_5_ori_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_6_ori_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_7_ori_20000.pdf"), width = 8, height = 6)
pdf(file = paste0(out_fp,"confusion_network_8_ori_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_network_1_newer_nobug_20000.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_retrain_single_dataset_process_mixture_strong.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_retrain_single_dataset_process_mixture_strong_better_testdata.pdf"), width = 8, height = 6)
# pdf(file = paste0(out_fp,"confusion_retrain_single_dataset_process_mixture_strong_known_param.pdf"), width = 8, height = 6)
print(plot_p)
# print(plot_2by2)
dev.off()
