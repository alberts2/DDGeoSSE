# # This R script takes complete trees (trees with extinct taxa) from julia simulation and drop extinct taxa from each tree.

  ###################
  # # LOAD PACKAGES
  ####################
  library(ape)
  
  drop_extinct_taxa = function(num_sim,batch_num){
    ###################
    # # SET FILE SYSTEM
    ####################
    fp     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/data/"
    in_fp  = paste0(fp,"sim",batch_num,"/complete_tree/") #input complete tree file from sim"batch_num"
    out_fp = paste0(fp,"sim",batch_num,"/extant_tree/") #output extant only tree for sim"batch_num"
    dat_fp = paste0(fp,"sim",batch_num,"/tip_data/") #chracter data for sim"batch_num"
    
    ###################
    # # READ TREE DATA
    ####################
    #
    not_exist = c()
    numsim = num_sim # total number of simulated trees 
    #
    for (i in 1:numsim){
      # sim_file = paste0(in_fp,"complete_tree_",i,".tre") # complete tree file name
      sim_file = paste0(in_fp,"sim.",i,".tre") # complete tree file name
      file.exists(sim_file) # check if tree i exists in file system
      #
      if (file.exists(sim_file)==TRUE){
        complete_tree = read.tree(file = sim_file) # read the complete tree file
        #
        # extant_tips   = read.csv(file = paste0(dat_fp,"tip_data_tree_",i,".dat.csv"),header = FALSE)[,1] # read the list of extant taxa
        extant_tips   = read.csv(file = paste0(dat_fp,"sim.",i,".dat.csv"),header = FALSE)[,1] # read the list of extant taxa
        # for (j in 1:length(extant_tips)){
        #   extant_tips[j] = paste0("sp.",extant_tips[j]) # rename with sp.
        # }
        extinct_tips = setdiff(complete_tree$tip.label,extant_tips) # get the list of extinct taxa
        #
        extant_tree = drop.tip(complete_tree,extinct_tips) # drop extinct taxa
        #
        write.tree(extant_tree,file=paste0(out_fp,"extant_tree_",i,".tre")) # save the reconstructed tree
        #
        print(paste0("extant tree ",i," done"))
      }
      else{
        not_exist = append(not_exist,i)
        print(paste0("complete tree ",i," does not exist"))
      }
    }
    # return(not_exist)
  }
  

