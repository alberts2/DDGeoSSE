# this script is to solve for n_i* given that we provide rate parameters and area size 
# (for another model variant that could give different n_i*)

library(nleqslv)

# Give set of parameters, find n_i_start for each region under unconstrained model

rho_w = 0.04 #base rate for speciation
rho_d = 0.002 #base rate for dispersal
rho_e = 0.001 #base rate for extinction

rho_w = exp(fit$par[1])
rho_d = exp(fit$par[2])
rho_e = exp(fit$par[3])

wD      = -2.5 #DD effect for speciation
dD_dest = -1.1 #DD effect for inbound dispersal
dD_src  = -2.2 #DD effect for outbound dispersal
eD      = 2 #DD effect for extinction

wD = fit$par[4]
dD_dest = fit$par[5]
dD_src = fit$par[6]
eD = fit$par[7]

# use log_n for the search so it does not fail if it searches negative values
equil_fun_log <- function(log_n) {
  n   = exp(log_n) #vector of ni* to find
  L   = log(n) + 1
  L_2 = log(n+1)+1
  K   = length(n)

  F = numeric(K)

  # balance equation for each region
  for (i in seq_len(K)) {
    incoming = sum(L[-i]^dD_src)

    F[i] = rho_w * L[i]^wD + rho_d * L_2[i]^dD_dest * incoming -rho_e * L[i]^eD
  }
  F
}


# fit <- nleqslv(
#   x = log(c(1, 10, 20)), #starting value for each ni_star
#   fn = equil_fun_log
# )

fit2 <- nleqslv(
  x = log(c(40, 20, 80)), #starting value for each ni_star
  fn = equil_fun_log
)

# with area size
rho_w = 0.04 #base rate for speciation
rho_d = 0.002 #base rate for dispersal
rho_e = 0.001 #base rate for extinction

wD      = -1.0 #DD effect for speciation
dD_dest = -1.0 #DD effect for inbound dispersal
dD_src  = 0 #DD effect for outbound dispersal
eD      = 1 #DD effect for extinction

equil_fun_log <- function(log_n) {
  n   = exp(log_n) #vector of ni* to find
  area_size = c(0.01,0.05,1) #relative area size 
  L   = log(n/area_size) + 1
  L_2 = log(n/area_size+1)+1
  K   = length(n)
  
  F = numeric(K)
  
  # balance equation for each region
  for (i in seq_len(K)) {
    incoming = sum(L[-i]^dD_src)
    
    F[i] = rho_w * L[i]^wD + rho_d * L_2[i]^dD_dest * incoming -rho_e * L[i]^eD
  }
  F
}

fit2 <- nleqslv(
  x = log(c(40, 20, 80)), #starting value for each ni_star
  fn = equil_fun_log
)
