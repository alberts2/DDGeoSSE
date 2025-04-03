files <- list.files(pattern = "^sim\\.[0-9]+\\.labels\\.csv$", full.names = FALSE)

indexes <- as.numeric(gsub("sim\\.([0-9]+)\\.labels\\.csv", "\\1", files))

time_vec = 0
num_vec = 0
for (i in 1:length(indexes)){
  dat=read.csv(paste0("sim.",indexes[i],".labels.csv"))
  if (dat$stop_condition == "numtaxa"){
    num_vec = num_vec+1
  }
  else {
    time_vec = time_vec+1
  }
}

prop_num = num_vec/(num_vec+time_vec)*100
prop_time = 100-prop_num