# this script gets the estimate for the same test datasets used on the paper (taken from 
# 400000 trees), make sure each model has equal representations in the test data, and within each model
# there's equal representations for trees with loguniform on base rates vs trees with uniform base rates
# then use these estimates to create the 8x3 confusion matrix 

# fp      = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/model_selection/Anolis/selection_sum_no_weight/"
fp      = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/deep_learning/server/DD-GeoSSE/5_reg/FINAL/mixture_dataset/model_selection/Anolis/selection_sum_weight_0.001/"
# ed_fp   = paste0(fp,"ed_1_anolis/estimate/")
# wd_fp   = paste0(fp,"wd_1_anolis/estimate/")
# dd_fp   = paste0(fp,"dd_1_anolis/estimate/")
# #
# ed_fp   = paste0(fp,"ed_2_anolis/estimate/")
# wd_fp   = paste0(fp,"wd_2_anolis/estimate/")
# dd_fp   = paste0(fp,"dd_2_anolis/estimate/")
# #
# ed_fp   = paste0(fp,"ed_3_anolis/estimate/")
# wd_fp   = paste0(fp,"wd_3_anolis/estimate/")
# dd_fp   = paste0(fp,"dd_3_anolis/estimate/")
# #
# ed_fp   = paste0(fp,"ed_4_anolis/estimate/")
# wd_fp   = paste0(fp,"wd_4_anolis/estimate/")
# dd_fp   = paste0(fp,"dd_4_anolis/estimate/")
# #
# ed_fp   = paste0(fp,"ed_5_anolis/estimate/")
# wd_fp   = paste0(fp,"wd_5_anolis/estimate/")
# dd_fp   = paste0(fp,"dd_5_anolis/estimate/")
# #
# ed_fp   = paste0(fp,"ed_6_anolis/estimate/")
# wd_fp   = paste0(fp,"wd_6_anolis/estimate/")
# dd_fp   = paste0(fp,"dd_6_anolis/estimate/")
# #
# ed_fp   = paste0(fp,"ed_7_anolis/estimate/")
# wd_fp   = paste0(fp,"wd_7_anolis/estimate/")
# dd_fp   = paste0(fp,"dd_7_anolis/estimate/")
#
ed_fp   = paste0(fp,"ed_8_anolis/estimate/")
wd_fp   = paste0(fp,"wd_8_anolis/estimate/")
dd_fp   = paste0(fp,"dd_8_anolis/estimate/")

# out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_1_ori_20000/"
# out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_2_ori_20000/"
# out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_3_ori_20000/"
# out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_4_ori_20000/"
# out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_5_ori_20000/"
# out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_6_ori_20000/"
# out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_7_ori_20000/"
out_fp  = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/data/validation_confusion/network_8_ori_20000/"

# all 20000 test samples for eD network
ed_all  = read.table(paste0(ed_fp,"sim.test_est.labels_cat.csv"),sep = ",",header = T)
ed_all_true  = read.table(paste0(ed_fp,"sim.test_true.labels_cat.csv"),sep = ",",header = T)

# all 20000 test samples for wD network
wd_all  = read.table(paste0(wd_fp,"sim.test_est.labels_cat.csv"),sep = ",",header = T)
wd_all_true  = read.table(paste0(wd_fp,"sim.test_true.labels_cat.csv"),sep = ",",header = T)

# all 20000 test samples for wD network
dd_all  = read.table(paste0(dd_fp,"sim.test_est.labels_cat.csv"),sep = ",",header = T)
dd_all_true  = read.table(paste0(dd_fp,"sim.test_true.labels_cat.csv"),sep = ",",header = T)

# divide dataset into test data that belong to each submodel 
wd_submodel_0 = wd_all[wd_all$idx > 100000 & wd_all$idx <= 150000,]
ed_submodel_0 = ed_all[ed_all$idx > 100000 & ed_all$idx <= 150000,]
dd_submodel_0 = dd_all[dd_all$idx > 100000 & dd_all$idx <= 150000,]
#
wd_submodel_0_true = wd_all_true[wd_all_true$idx > 100000 & wd_all_true$idx <= 150000,]
ed_submodel_0_true = ed_all_true[ed_all_true$idx > 100000 & ed_all_true$idx <= 150000,]
dd_submodel_0_true = dd_all_true[dd_all_true$idx > 100000 & dd_all_true$idx <= 150000,]

write.table(wd_submodel_0,paste0(out_fp,"wd_1/estimate_0/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_0_true,paste0(out_fp,"wd_1/estimate_0/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_0,paste0(out_fp,"ed_1/estimate_0/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_0_true,paste0(out_fp,"ed_1/estimate_0/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_0,paste0(out_fp,"dd_1/estimate_0/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_0_true,paste0(out_fp,"dd_1/estimate_0/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

wd_submodel_1 = wd_all[wd_all$idx > 0 & wd_all$idx <= 50000,]
ed_submodel_1 = ed_all[ed_all$idx > 0 & ed_all$idx <= 50000,]
dd_submodel_1 = dd_all[dd_all$idx > 0 & dd_all$idx <= 50000,]
#
wd_submodel_1_true = wd_all_true[wd_all_true$idx > 0 & wd_all_true$idx <= 50000,]
ed_submodel_1_true = ed_all_true[ed_all_true$idx > 0 & ed_all_true$idx <= 50000,]
dd_submodel_1_true = dd_all_true[dd_all_true$idx > 0 & dd_all_true$idx <= 50000,]

write.table(wd_submodel_1,paste0(out_fp,"wd_1/estimate_1/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_1_true,paste0(out_fp,"wd_1/estimate_1/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_1,paste0(out_fp,"ed_1/estimate_1/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_1_true,paste0(out_fp,"ed_1/estimate_1/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_1,paste0(out_fp,"dd_1/estimate_1/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_1_true,paste0(out_fp,"dd_1/estimate_1/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

wd_submodel_2 = wd_all[wd_all$idx > 50000 & wd_all$idx <= 100000,]
ed_submodel_2 = ed_all[ed_all$idx > 50000 & ed_all$idx <= 100000,]
dd_submodel_2 = dd_all[dd_all$idx > 50000 & dd_all$idx <= 100000,]
#
wd_submodel_2_true = wd_all_true[wd_all_true$idx > 50000 & wd_all_true$idx <= 100000,]
ed_submodel_2_true = ed_all_true[ed_all_true$idx > 50000 & ed_all_true$idx <= 100000,]
dd_submodel_2_true = dd_all_true[dd_all_true$idx > 50000 & dd_all_true$idx <= 100000,]

write.table(wd_submodel_2,paste0(out_fp,"wd_1/estimate_2/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_2_true,paste0(out_fp,"wd_1/estimate_2/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_2,paste0(out_fp,"ed_1/estimate_2/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_2_true,paste0(out_fp,"ed_1/estimate_2/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_2,paste0(out_fp,"dd_1/estimate_2/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_2_true,paste0(out_fp,"dd_1/estimate_2/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

wd_submodel_3 = wd_all[wd_all$idx > 150000 & wd_all$idx <= 200000,]
ed_submodel_3 = ed_all[ed_all$idx > 150000 & ed_all$idx <= 200000,]
dd_submodel_3 = dd_all[dd_all$idx > 150000 & dd_all$idx <= 200000,]
#
wd_submodel_3_true = wd_all_true[wd_all_true$idx > 150000 & wd_all_true$idx <= 200000,]
ed_submodel_3_true = ed_all_true[ed_all_true$idx > 150000 & ed_all_true$idx <= 200000,]
dd_submodel_3_true = dd_all_true[dd_all_true$idx > 150000 & dd_all_true$idx <= 200000,]
#
write.table(wd_submodel_3,paste0(out_fp,"wd_1/estimate_3/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_3_true,paste0(out_fp,"wd_1/estimate_3/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_3,paste0(out_fp,"ed_1/estimate_3/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_3_true,paste0(out_fp,"ed_1/estimate_3/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_3,paste0(out_fp,"dd_1/estimate_3/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_3_true,paste0(out_fp,"dd_1/estimate_3/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

wd_submodel_4 = wd_all[wd_all$idx > 200000 & wd_all$idx <= 250000,]
ed_submodel_4 = ed_all[ed_all$idx > 200000 & ed_all$idx <= 250000,]
dd_submodel_4 = dd_all[dd_all$idx > 200000 & dd_all$idx <= 250000,]
#
wd_submodel_4_true = wd_all_true[wd_all_true$idx > 200000 & wd_all_true$idx <= 250000,]
ed_submodel_4_true = ed_all_true[ed_all_true$idx > 200000 & ed_all_true$idx <= 250000,]
dd_submodel_4_true = dd_all_true[dd_all_true$idx > 200000 & dd_all_true$idx <= 250000,]

write.table(wd_submodel_4,paste0(out_fp,"wd_1/estimate_4/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_4_true,paste0(out_fp,"wd_1/estimate_4/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_4,paste0(out_fp,"ed_1/estimate_4/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_4_true,paste0(out_fp,"ed_1/estimate_4/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_4,paste0(out_fp,"dd_1/estimate_4/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_4_true,paste0(out_fp,"dd_1/estimate_4/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

wd_submodel_5 = wd_all[wd_all$idx > 250000 & wd_all$idx <= 300000,]
ed_submodel_5 = ed_all[ed_all$idx > 250000 & ed_all$idx <= 300000,]
dd_submodel_5 = dd_all[dd_all$idx > 250000 & dd_all$idx <= 300000,]
#
wd_submodel_5_true = wd_all_true[wd_all_true$idx > 250000 & wd_all_true$idx <= 300000,]
ed_submodel_5_true = ed_all_true[ed_all_true$idx > 250000 & ed_all_true$idx <= 300000,]
dd_submodel_5_true = dd_all_true[dd_all_true$idx > 250000 & dd_all_true$idx <= 300000,]

write.table(wd_submodel_5,paste0(out_fp,"wd_1/estimate_5/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_5_true,paste0(out_fp,"wd_1/estimate_5/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_5,paste0(out_fp,"ed_1/estimate_5/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_5_true,paste0(out_fp,"ed_1/estimate_5/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_5,paste0(out_fp,"dd_1/estimate_5/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_5_true,paste0(out_fp,"dd_1/estimate_5/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

wd_submodel_6 = wd_all[wd_all$idx > 300000 & wd_all$idx <= 350000,]
ed_submodel_6 = ed_all[ed_all$idx > 300000 & ed_all$idx <= 350000,]
dd_submodel_6 = dd_all[dd_all$idx > 300000 & dd_all$idx <= 350000,]
#
wd_submodel_6_true = wd_all_true[wd_all_true$idx > 300000 & wd_all_true$idx <= 350000,]
ed_submodel_6_true = ed_all_true[ed_all_true$idx > 300000 & ed_all_true$idx <= 350000,]
dd_submodel_6_true = dd_all_true[dd_all_true$idx > 300000 & dd_all_true$idx <= 350000,]

write.table(wd_submodel_6,paste0(out_fp,"wd_1/estimate_6/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_6_true,paste0(out_fp,"wd_1/estimate_6/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_6,paste0(out_fp,"ed_1/estimate_6/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_6_true,paste0(out_fp,"ed_1/estimate_6/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_6,paste0(out_fp,"dd_1/estimate_6/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_6_true,paste0(out_fp,"dd_1/estimate_6/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

wd_submodel_7 = wd_all[wd_all$idx > 350000 & wd_all$idx <= 400000,]
ed_submodel_7 = ed_all[ed_all$idx > 350000 & ed_all$idx <= 400000,]
dd_submodel_7 = dd_all[dd_all$idx > 350000 & dd_all$idx <= 400000,]
#
wd_submodel_7_true = wd_all_true[wd_all_true$idx > 350000 & wd_all_true$idx <= 400000,]
ed_submodel_7_true = ed_all_true[ed_all_true$idx > 350000 & ed_all_true$idx <= 400000,]
dd_submodel_7_true = dd_all_true[dd_all_true$idx > 350000 & dd_all_true$idx <= 400000,]

write.table(wd_submodel_7,paste0(out_fp,"wd_1/estimate_7/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(wd_submodel_7_true,paste0(out_fp,"wd_1/estimate_7/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(ed_submodel_7,paste0(out_fp,"ed_1/estimate_7/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(ed_submodel_7_true,paste0(out_fp,"ed_1/estimate_7/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
#
write.table(dd_submodel_7,paste0(out_fp,"dd_1/estimate_7/sim.test_est.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(dd_submodel_7_true,paste0(out_fp,"dd_1/estimate_7/sim.test_true.labels_cat.csv"),sep=",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
