# # This R script takes complete trees (trees with extinct taxa) from julia simulation and drop extinct taxa from each tree.

  ###################
  # # LOAD PACKAGES
  ####################
  library(ape)
  
  ###################
  # # SET FILE SYSTEM
  ####################
  fp = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/code/jl/example/"
  
  ###################
  # # READ TREE DATA
  ####################
  #
  numsim = 10 # total number of simulated trees 
  #
  for (i in 1:numsim){
    sim_file = paste0(fp,"complete_tree_",i,".tre") # tree file name
    file.exists(sim_file) # check if tree i exists in file system
    #
    if (file.exists(sim_file)==TRUE){
      complete_tree = read.tree(file = sim_file) # read the tree file
      #
      extant_tips   = read.csv(file = paste0(fp,"tip_data_tree_",i),header = FALSE)[,1] # read the list of extant taxa
      for (j in 1:length(extant_tips)){
        extant_tips[j] = paste0("sp.",extant_tips[j]) # rename with sp.
      }
      extinct_tips = setdiff(complete_tree$tip.label,extant_tips) # get the list of extinct taxa
      #
      extant_tree = drop.tip(complete_tree,extinct_tips) # drop extinct taxa
      #
      write.tree(extant_tree,file=paste0(fp,"extant_tree_",i,".tre")) # save the reconstructed tree
    }
  }


