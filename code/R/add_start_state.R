# This code adds starting state for each node -- useful for ancestral state recons with phyddle
# where we encode end state of parent and start state of both childrens as a joint state

in_fp = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/code/jl/simulate/"
setwd(in_fp)

num_sim = 10

for (i in 1:num_sim){
  file_paths = paste0("cladogenetic event history_",i-1,".csv")
  if (length(file_paths[file.exists(file_paths)]) == 0){
    print(paste0("simulation ", i-1, " missing"))
  } else {
    clado_dat   = read.table(file_paths,header = TRUE, sep = ",", stringsAsFactors = FALSE)
    anagen_dat  = read.table(paste0("anagenetic_changes_tree_",i-1,".csv"),header = TRUE, sep = ",", stringsAsFactors = FALSE)
    #
    num_repli = rep(NA,nrow(clado_dat))
    full_history = data.frame(parent_idx = num_repli, child_1_idx =  num_repli, child_2_idx =  num_repli,
                              parent_start_state = num_repli, parent_end_state = num_repli,
                              child_1_start_state = num_repli, child_1_end_state = num_repli,
                              child_2_start_state = num_repli, child_2_end_state = num_repli,
                              extinct = num_repli, split = num_repli,
                              start_time = num_repli, end_time = num_repli, branch_length = num_repli)
    #
    for (j in 1:nrow(clado_dat)){ #loop over all nodes
      parent_idx = clado_dat$idx[j] #get the parent index
      child1_idx = clado_dat$child_1_idx[j] #get child 1 index
      child2_idx = clado_dat$child_2_idx[j] #get child 2 index
      #
      full_history$parent_idx[j]       = parent_idx #append parent index
      full_history$child_1_idx[j]      = child1_idx #append child 1 index
      full_history$child_2_idx[j]      = child2_idx #append child 2 index
      full_history$extinct[j]          = clado_dat$extinct_status[j]
      full_history$split[j]            = clado_dat$split_status[j]
      full_history$start_time[j]       = clado_dat$time_start[j]
      full_history$end_time[j]         = clado_dat$time_end[j]
      full_history$branch_length[j]    = clado_dat$branch_length[j]
      # get start and end state of parent idx
      which_parent_idx = which(anagen_dat$idx==parent_idx)
      if (length(which_parent_idx) == 0){
        full_history$parent_start_state[j] = clado_dat$final_state[parent_idx]
        full_history$parent_end_state[j]   = full_history$parent_start_state[j] 
      } else {
        full_history$parent_start_state[j] = anagen_dat$from_state[which_parent_idx[1]]
        full_history$parent_end_state[j]   = clado_dat$final_state[parent_idx]
      }
      # get start and end state of child 1 idx
      which_child1_idx = which(anagen_dat$idx==child1_idx)
      if (length(which_child1_idx) == 0){
        full_history$child_1_start_state[j] = clado_dat$final_state[child1_idx]
        full_history$child_1_end_state[j] = full_history$child_1_start_state[j]
      } else {
        full_history$child_1_start_state[j] = anagen_dat$from_state[which_child1_idx[1]]
        full_history$child_1_end_state[j]   = clado_dat$final_state[child1_idx]
      }
      # get start and end state of child 2 idx
      which_child2_idx = which(anagen_dat$idx==child2_idx)
      if (length(which_child2_idx) == 0){
        full_history$child_2_start_state[j] = clado_dat$final_state[child2_idx]
        full_history$child_2_end_state[j] = full_history$child_2_start_state[j]
      } else {
        full_history$child_2_start_state[j] = anagen_dat$from_state[which_child2_idx[1]]
        full_history$child_2_end_state[j]   = clado_dat$final_state[child2_idx]
      }
    }
  }
  write.table(full_history, paste0("full_history_", i - 1, ".csv"), sep = ",",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}
