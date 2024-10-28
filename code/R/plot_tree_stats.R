# This script returns summary statistics across set of trees such as 
# -. mean number of extant species
# -. mean number of extant + extinct species 
# -. mean number of event per type (w, local extinction, global extinction, d, and, b)
# and plot barchart displaying the proportion for each event type for each tree. 

library(ape)
library(ggplot2)

plot_tree_stats = function(num_sim,batch_num){
  ###################
  # # SET FILE SYSTEM
  ###################
  fp     = "/Users/albertsoewongsono/Documents/Code Testing/density_geosse/"
  tree_dat_fp     = paste0(fp,"data/sim",batch_num,"/complete_tree/") # complete tree data
  extant_dat_fp   = paste0(fp,"data/sim",batch_num,"/extant_tree/") # extant-only tree data 
  anagen_dat_fp   = paste0(fp,"data/sim",batch_num,"/anagen_changes/") # anagenetic event data 
  cladogen_dat_fp = paste0(fp,"data/sim",batch_num,"/cladogen_history/") # cladogenetic event data
  tipstate_dat_fp = paste0(fp,"data/sim",batch_num,"/tip_data/") # tip state data
  
  out_fp = paste0(fp,"plot/sim",batch_num,"/") #output 
  
  ###################
  # # INITIALIZATION
  ####################
  num_extant_tips_vec = c() # vector containing number of extant taxa across num_sim simulations 
  num_total_tips_vec  = c() # vector containing number of total taxa (including extinct) across num_sim simulations 
  #
  num_dispersal_vec   = c() # vector containing number of dispersal events across num_sim simulations
  num_within_vec      = c() # vector containing number of within-region speciation events across num_sim simulations
  num_between_vec     = c() # vector containing number of between-region speciation events across num_sim simulations
  num_localext_vec    = c() # vector containing number of local extinction events across num_sim simulations
  num_globalext_vec   = c() # vector containing number of global extinction events across num_sim simulations
  #
  individual_stats    = c() # storing each event count for each tree e.g. c(tree 1 stats,tree 2 stats)
  #
  num_sim    = num_sim # number of simulated trees
  
  for (i in 1:num_sim){
    ###################
    # # READ TREE DATA
    ####################
    complete_tree = read.tree(paste0(tree_dat_fp,"complete_tree_",i,".tre")) # read complete tree i
    extant_tree   = read.tree(paste0(extant_dat_fp,"extant_tree_",i,".tre")) # read extant-only tree i
    ###################
    # # READ EVENT DATA
    ####################
    anagen_df   = read.table(paste0(anagen_dat_fp,"anagenetic_changes_tree_",i,".csv"),sep = ",",header = T) # read anagenetic event history for tree i
    cladogen_df = read.table(paste0(cladogen_dat_fp,"cladogenetic event history_",i,".csv"),sep = ",",header = T) # read cladogenetic event history for tree i
    #
    # append to total tips and extant tips for tree i 
    num_extant_tips_vec = append(num_extant_tips_vec,length(extant_tree$tip.label))
    num_total_tips_vec  = append(num_total_tips_vec,length(complete_tree$tip.label))
    # count each event type for tree i
    count_dispersal     = nrow(subset(anagen_df,from_state < to_state))
    count_localext      = nrow(subset(anagen_df,from_state > to_state))
    count_within        = nrow(subset(cladogen_df,event_type == "w"))
    count_between       = nrow(subset(cladogen_df,event_type == "b"))
    count_globalext     = nrow(subset(cladogen_df,event_type == "e"))
    # append total number of each event type for tree i
    num_dispersal_vec   = append(num_dispersal_vec,count_dispersal)
    num_localext_vec    = append(num_localext_vec,count_localext)
    num_within_vec      = append(num_within_vec,count_within)
    num_between_vec     = append(num_between_vec,count_between)
    num_globalext_vec   = append(num_globalext_vec,count_globalext)
    # append each event type count for tree i
    individual_stats    = append(individual_stats,c(count_within,count_between,count_dispersal,count_localext,count_globalext))
    # print some message about simulation i data
    print(paste0("simulation ",i," stats are done"))
  }
  #
  ###################
  # # OUTPUT
  ####################
  # dataframe for summarizing mean statistic for the event count data across num_sim trees 
  summary_df = data.frame("total_trees" = num_sim,"mean_taxa_count" = mean(num_total_tips_vec),"mean_extant_taxa_count" = mean(num_extant_tips_vec),
                          "mean_within_region_count" = mean(num_within_vec),"mean_between_region_count" = mean(num_between_vec),
                          "mean_dispersal_count" = mean(num_dispersal_vec),"mean_local_ext_count" = mean(num_localext_vec),
                          "mean_global_ext_count" = mean(num_globalext_vec))
  #
  # dataframe for stacked barchart displaying proportion of each event type for each tree 
  plot_df = data.frame(tree_id = rep(seq(1,num_sim),each = 5),event_types = rep(c("within-region","between-region","dispersal","local ext","global ext"),num_sim),
                       event_counts = individual_stats )
  ###################
  # # PLOTS (as stacked barplot for number of different events for each tree?)
  ## also run again simtest 11 on server for simpler case, see if it still works 
  ####################
  #
  p = ggplot(plot_df, aes(fill=event_types, y=event_counts, x=tree_id)) + 
      geom_bar(position="fill", stat="identity") + 
      scale_fill_manual(name = "Event type",values=c('within-region'='#3E9AB4','between-region'='#DA9556','dispersal'='#B2DE81','local ext'='#BFC2C1','global ext'='#E03A42')) +
      labs(title = "Proportion of event types across individual trees", x = "Tree ID", y = "Event proportion") +
      theme(plot.title = element_text(hjust = 0.5))
  #
  plot(p)
  ###################
  # # SAVE OUTPUT
  ####################
  write.csv(summary_df,file = paste0(out_fp,"summary_statistics.csv"),row.names = F) # save the summary statistics
  ggsave("event_type_proportion.jpg",plot = p, path = out_fp,width = 8, height = 6, units = "in") # save the stacked barchart
}
  