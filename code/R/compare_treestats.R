# This script compares tree balance and branch length statistics between julia and diversitree simulations 
# for standard GeoSSE model

library(treestats)
library(ape)

##### FILE SYSTEM #####
#######################
fp       = "/Users/albertsoewongsono/Documents/Code\ Testing/density_geosse/"
in_fp    = paste0(fp,"data/sim_geosse/diversitree/") 
in_jl_fp = paste0(fp,"data/sim_geosse/julia/")
# in_jl_fp = paste0(fp,"data/sim_geosse/julia_false/")
out_fp   = paste0(fp,"plot/sim_geosse/")
# out_fp   = paste0(fp,"plot/sim_geosse_false/")

##### INPUT #####
#################

numsim = 1000

# Read diversitree input 
diversitree_dat = read.table(file = paste0(in_fp,"diversitree_dat.csv"),header = T,sep = ",") # diversitree's topology data 
diversitree_length_dat = read.table(file = paste0(in_fp,"diversitree_length_dat.csv"),header = T,sep = ",") # diversitree's branch length data 

# Read Julia input 
balance_jl_df = data.frame(colless = NA, colless_corr = NA, colless_quad = NA,
                           sackin = NA, beta = NA, stairs = NA, stairs2 = NA,
                           blum_norm = NA, rogers_norm = NA, pitchforks_norm = NA)


length_jl_df = data.frame(crown_age = NA, gamma = NA, pigot_rho = NA)

for (i in 1:1000){
  phy_dat = read.tree(file = paste0(in_jl_fp,"extant_tree/","extant_tree_",i,".tre"))
  # remove the stem branch to match with diversitree simulation
  phy_dat$root.edge = 0
  #
  if (phy_dat$Nnode > 2){ #ignore extant-only julia trees with <= 2 tips since can't compute beta with only 2 taxa
    balance_jl_df[i,] = c(colless(phy_dat),colless_corr(phy_dat),colless_quad(phy_dat),
                          sackin(phy_dat),beta_statistic(phy_dat),stairs(phy_dat),stairs2(phy_dat),
                          blum(phy_dat,normalization = T),rogers(phy_dat,normalization = "tips"),
                          pitchforks(phy_dat,normalization = "tips"))
    #
    length_jl_df[i,]  = c(crown_age(phy_dat),treestats::gamma_statistic(phy_dat),pigot_rho(phy_dat))
  }
  balance_jl_df = na.omit(balance_jl_df) # remove rows with na values
  length_jl_df  = na.omit(length_jl_df)  # remove rows with na values
}

diversitree_dat = diversitree_dat[-(1:(numsim-nrow(balance_jl_df))),] #adjust numsim accordingly to match num rows from julia dataframe
diversitree_length_dat = diversitree_length_dat[-(1:(numsim-nrow(length_jl_df))),] #adjust numsim accordingly to match num rows from julia dataframe


##### SUMMARY ###
#################

# colless stat
colless_dat = data.frame(colless_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                        value = c(diversitree_dat$colless,balance_jl_df$colless))

# colless_corr stat
colless_corr_dat = data.frame(colless_corr_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                              value = c(diversitree_dat$colless_corr,balance_jl_df$colless_corr))

# colless_quad stat
colless_quad_dat = data.frame(colless_quad_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                              value = c(diversitree_dat$colless_quad,balance_jl_df$colless_quad))

# sackin stat
sackin_dat = data.frame(sackin_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                              value = c(diversitree_dat$sackin,balance_jl_df$sackin))

# beta stat
beta_dat = data.frame(beta_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                        value = c(diversitree_dat$beta,balance_jl_df$beta))

# stairs stat
stairs_dat = data.frame(stairs_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                      value = c(diversitree_dat$stairs,balance_jl_df$stairs))

# stairs2 stat
stairs2_dat = data.frame(stairs2_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                        value = c(diversitree_dat$stairs2,balance_jl_df$stairs2))


# blum stat
blum_dat = data.frame(blum_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                      value = c(diversitree_dat$blum_norm,balance_jl_df$blum_norm))

# rogers stat
rogers_dat = data.frame(rogers_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                        value = c(diversitree_dat$rogers_norm,balance_jl_df$rogers_norm))

# pitchforks stat
pitchforks_dat = data.frame(pitchforks_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                            value = c(diversitree_dat$pitchforks_norm,balance_jl_df$pitchforks_norm))

# crown age stat
crown_dat = data.frame(crown_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                         value = c(diversitree_length_dat$crown_age,length_jl_df$crown_age))

# gamma stat
gamma_dat = data.frame(gamma_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                       value = c(diversitree_length_dat$gamma,length_jl_df$gamma))

# pigot_rho stat
rho_dat = data.frame(rho_type=sort(rep(c("julia","diversitree"),nrow(diversitree_dat))), 
                       value = c(diversitree_length_dat$pigot_rho,length_jl_df$pigot_rho))


##### BALANCE PLOT ###
#################

plot_colless_corr = ggplot(colless_corr_dat, aes(x=colless_corr_type, y=value, fill=colless_corr_type)) +
                    geom_boxplot() +
                    ggtitle('colless_corr')  +  scale_fill_discrete(name = "colless_corr") + 
                    labs(x = "source") + ylim(c(0,0.04)) 

plot_colless = ggplot(colless_dat, aes(x=colless_type, y=value, fill=colless_type)) +
               geom_boxplot() +
               ggtitle('colless')  + scale_fill_discrete(name = "colless") + 
               labs(x = "source") 

plot_colless_quad = ggplot(colless_quad_dat, aes(x=colless_quad_type, y=value, fill=colless_quad_type)) +
                    geom_boxplot() +
                    ggtitle('colless_quad')  + scale_fill_discrete(name = "colless_quad") + 
                    labs(x = "source") + ylim(c(0,1e+07))

plot_sackin = ggplot(sackin_dat, aes(x=sackin_type, y=value, fill=sackin_type)) +
              geom_boxplot() +
              ggtitle('sackin') + scale_fill_discrete(name = "sackin") + 
              labs(x = "source") 

plot_beta = ggplot(beta_dat, aes(x=beta_type, y=value, fill=beta_type)) +
            geom_boxplot() +
            ggtitle('beta') + scale_fill_discrete(name = "beta") + 
            labs(x = "source") + ylim(c(0,1))

plot_stairs = ggplot(stairs_dat, aes(x=stairs_type, y=value, fill=stairs_type)) +
              geom_boxplot() +
              ggtitle('stairs') + scale_fill_discrete(name = "stairs") + 
              labs(x = "source") + ylim(c(0.4,0.8))

plot_stairs2 = ggplot(stairs2_dat, aes(x=stairs2_type, y=value, fill=stairs2_type)) +
               geom_boxplot() + ggtitle('stairs2') +
               scale_fill_discrete(name = "stairs2") + 
               labs(x = "source") + ylim(c(0.5,0.7))

plot_blum = ggplot(blum_dat, aes(x=blum_type, y=value, fill=blum_type)) +
            geom_boxplot() +
            ggtitle('normalized blum') + scale_fill_discrete(name = "blum") + 
            labs(x = "source") 

plot_rogers = ggplot(rogers_dat, aes(x=rogers_type, y=value, fill=rogers_type)) +
              geom_boxplot() +
              ggtitle('normalized rogers') + scale_fill_discrete(name = "rogers") + 
              labs(x = "source") 

plot_pitchforks = ggplot(pitchforks_dat, aes(x=pitchforks_type, y=value, fill=pitchforks_type)) +
                  geom_boxplot() +
                  ggtitle('normalized pitchforks') + scale_fill_discrete(name = "pitchforks") + 
                  labs(x = "source") 

##### BRANCH LENGTH PLOT ###
#################

plot_crown = ggplot(crown_dat, aes(x=crown_type, y=value, fill=crown_type)) +
  geom_boxplot() +
  ggtitle('crown_age')  + scale_fill_discrete(name = "crown age") + 
  labs(x = "source") 

plot_gamma = ggplot(gamma_dat, aes(x=gamma_type, y=value, fill=gamma_type)) +
  geom_boxplot() +
  ggtitle('gamma')  + scale_fill_discrete(name = "gamma") + 
  labs(x = "source") 

plot_rho = ggplot(rho_dat, aes(x=rho_type, y=value, fill=rho_type)) +
  geom_boxplot() +
  ggtitle('pigot_rho')  + scale_fill_discrete(name = "pigot rho") + 
  labs(x = "source") 

##### SAVE BALANCE OUTPUT ###
#####################

pdf(file = paste0(out_fp,"comparison_geosse.pdf"), height = 6, width = 8)

print(plot_colless_corr)
print(plot_colless)
print(plot_colless_quad)
print(plot_sackin)
print(plot_beta)
print(plot_stairs)
print(plot_stairs2)
print(plot_blum)
print(plot_rogers)
print(plot_pitchforks)

dev.off()

##### SAVE BRANCH LENGTH OUTPUT ###
#####################

pdf(file = paste0(out_fp,"comparison_branch_geosse.pdf"), height = 6, width = 8)

print(plot_crown)
print(plot_gamma)
print(plot_rho)

dev.off()

##### T-test for sample means ###
#################################
p_value_vec = c(t.test(diversitree_dat$colless,balance_jl_df$colless)$p.value,
                t.test(diversitree_dat$colless_corr,balance_jl_df$colless_corr)$p.value,
                t.test(diversitree_dat$colless_quad,balance_jl_df$colless_quad)$p.value,
                t.test(diversitree_dat$sackin,balance_jl_df$sackin)$p.value,
                t.test(diversitree_dat$beta,balance_jl_df$beta)$p.value,
                t.test(diversitree_dat$stairs,balance_jl_df$stairs)$p.value,
                t.test(diversitree_dat$stairs2,balance_jl_df$stairs2)$p.value,
                t.test(diversitree_dat$blum_norm,balance_jl_df$blum_norm)$p.value,
                t.test(diversitree_dat$rogers_norm,balance_jl_df$rogers_norm)$p.value,
                t.test(diversitree_dat$pitchforks_norm,balance_jl_df$pitchforks_norm)$p.value)
