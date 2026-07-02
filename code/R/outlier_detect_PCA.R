# Identify outlier for aux data using LoF criterion on PCA space

# Load package 
library(rhdf5) # read rhdf5 file
library(reticulate) # to run python on R

# call venv for importing python's sklearn library
use_python("~/r-reticulate-sklearn/bin/python", required = TRUE)
sklearn_neighbors = import("sklearn.neighbors")

# File setting 
file_fn     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/data/sensitivity_test/"
# submodel  = paste0(file_fn,"submodel_1_sensitivity/format_weight/")
# submodel  = paste0(file_fn,"submodel_2_sensitivity/format_weight/")
# submodel  = paste0(file_fn,"submodel_3_sensitivity/format_weight/")
# submodel  = paste0(file_fn,"submodel_4_sensitivity/format_weight/")
# submodel  = paste0(file_fn,"submodel_5_sensitivity/format_weight/")
# submodel  = paste0(file_fn,"submodel_6_sensitivity/format_weight/")
submodel  = paste0(file_fn,"submodel_7_sensitivity/format_weight/")
setwd(submodel)

# names = h5ls("sim.empirical.hdf5")

# empirical data with relative area size
emp_df = as.data.frame(h5read("sim.empirical.hdf5", "/aux_data"))

# empirical feature name 
emp_name = h5read("sim.empirical.hdf5", "/aux_data_names")

rownames(emp_df) = emp_name

# train data without relative area size 
train_df = as.data.frame(h5read("sim.train.hdf5", "/aux_data"))

# train feature name
train_name = h5read("sim.train.hdf5", "/aux_data_names")

rownames(train_df) = train_name

# transpose feature as columns instead
emp_df   = as.data.frame(t(emp_df))
train_df = as.data.frame(t(train_df))

# make sure both emp and train data have the same features 
emp_df = emp_df[, colnames(train_df), drop = FALSE]

# LOF outlier detection
num_comp     = 4  # num of PCs
k            = 20 # The kth-distance to be used to calculate LOFs
# lof_quantile = 0.99995 
lof_quantile = 0.95

# remove any invalid  training features
train_mean0 = colMeans(train_df) 
train_sd0   = apply(train_df, 2, sd)

keep = is.finite(train_mean0) &
       is.finite(train_sd0) &
       train_sd0 > 0

train_df2 = train_df[, keep, drop = FALSE] 
emp_df2   = emp_df[, keep, drop = FALSE]

# standardize using training mean/sd
train_mean = colMeans(train_df2)  # feature mean from training data
train_sd   = apply(train_df2, 2, sd) # feature sd from training data

train_z = sweep(train_df2, 2, train_mean, "-")
train_z = sweep(train_z, 2, train_sd, "/")

# standardize empirical data using training mean and sd
emp_z = sweep(emp_df2, 2, train_mean, "-")
emp_z = sweep(emp_z, 2, train_sd, "/")

# PCA fitted on training data only
num_comp = min(num_comp, ncol(train_z)) # number of PCs

pca_fit = prcomp(train_z, center = FALSE, scale. = FALSE) # fit PCA on training

train_pca = pca_fit$x[, seq_len(num_comp), drop = FALSE] # get the train PC coordinates 

# project the empirical data into PCA space 
emp_pca = predict(pca_fit, newdata = emp_z)[, seq_len(num_comp), drop = FALSE]

# Use LoF outlier detection from sklearn python library
lof = sklearn_neighbors$LocalOutlierFactor(
        n_neighbors = as.integer(k), # scikit imports the number as a float, need to convert
        novelty = TRUE
)

lof$fit(as.matrix(train_pca)) # fit LoF model using train PC coordinates

lof_est      = lof$predict(as.matrix(emp_pca)) # classify each emp sample based on fitted lof model
lof_decision = lof$decision_function(as.matrix(emp_pca))

lof_summary = data.frame(
                sample         = rownames(emp_df2),
                prediction     = as.integer(lof_est),
                decision_score = as.numeric(lof_decision),
                outlier        = as.integer(lof_est) == -1
              )

prop_outlier = nrow(lof_summary[lof_summary$outlier=="TRUE",])/250*100
