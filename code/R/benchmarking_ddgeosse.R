# this script plots the benchmarking time across 1000 simulations with 
# 3 to 7 regions conditioned on number of extant taxa equals to 150. 
# In each simulation, we start with species in region A. 

library(ggplot2)

fp              = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
dat_fp          = paste0(fp,"data/benchmarking/")
out_fp          = paste0(fp,"plot/benchmarking/")

# load the simulated data

dat_3reg = read.table(paste0(dat_fp,"simulate_3reg/sim.benchmark_3reg.csv"),sep = ",",header = T)
dat_4reg = read.table(paste0(dat_fp,"simulate_4reg/sim.benchmark_4reg.csv"),sep = ",",header = T)
dat_5reg = read.table(paste0(dat_fp,"simulate_5reg/sim.benchmark_5reg.csv"),sep = ",",header = T)
dat_6reg = read.table(paste0(dat_fp,"simulate_6reg/sim.benchmark_6reg.csv"),sep = ",",header = T)
dat_7reg = read.table(paste0(dat_fp,"simulate_7reg/sim.benchmark_7reg.csv"),sep = ",",header = T)

time_3reg = dat_3reg$elapsed_sec
time_4reg = dat_4reg$elapsed_sec
time_5reg = dat_5reg$elapsed_sec
time_6reg = dat_6reg$elapsed_sec
time_7reg = dat_7reg$elapsed_sec

comb_dat = data.frame(num_reg = c(3,4,5,6,7),
                      median_time = c(median(time_3reg),
                                      median(time_4reg),
                                      median(time_5reg),
                                      median(time_6reg),
                                      median(time_7reg)))

comb_dat$region <- factor(
  comb_dat$num_reg,
  levels = c(3,4,5,6,7),
  labels = c("3 regions", "4 regions", "5 regions","6 regions","7 regions")
)

p = ggplot(comb_dat, aes(x = num_reg, y = median_time)) +
  geom_line(group = 1) +
  geom_point(aes(color = region), size = 5) +
  scale_x_continuous(breaks = c(3,4,5,6,7)) +
  theme_minimal() +
  labs(
    x = "Number of Regions",
    y = "Median Elapsed Time (in seconds)",
    color = NULL
  ) + 
  theme(axis.text.x = element_text(size = 24),   # x-axis tick labels
        axis.text.y = element_text(size = 24),   # y-axis tick labels
        axis.title.x = element_text(size = 24),  # x-axis title
        axis.title.y = element_text(size = 24),  # y-axis title
        legend.text = element_text(size = 16),   # legend item labels
        legend.title = element_text(size = 16)) 

# region 6 and 7 have no carrying capacity? check again


pdf(file = paste0(out_fp,"benchmarking_ddgeosse.pdf"), width = 8, height = 6)
print(p)
dev.off()

