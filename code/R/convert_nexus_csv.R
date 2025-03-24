base_path      = Sys.getenv("HOME") # get to user's base directory
fp             = file.path(base_path,"Documents","Code Testing","density_geosse","data","empiric_viburnum","viburnum_dat_edit")

setwd(fp)

nex_dat = read.nexus.data("oreino5.1.nex")

viburnum_dat = data.frame("species" = NA, "reg1" = NA, "reg2" = NA,
                          "reg3" = NA, "reg4" = NA, "reg5" = NA)

for (i in 1:length(nex_dat)){
  row_dat = c()
  row_dat = append(row_dat,names(nex_dat[i]))
  row_dat = append(row_dat,as.numeric(nex_dat[[i]]))
  #
  viburnum_dat = rbind(viburnum_dat,row_dat)
}
viburnum_dat = viburnum_dat[-1,]

write.csv(viburnum_dat,"oreino5.1.dat.csv",row.names = F)