  # Given target equilibrium diversity (unequal) across each region, find shared set of parameters
  # by optimizing balance equation across regions 
  
  ## Using the current formulation without area size (unique to each region)
  ## This causes exchangebility problem in which we cannot identify which region with which 
  ## equilibrium diversity, in case this parameter (equilibrium) in not equal across all regions 
  ## this is because n_i, although it gives region-dep rates, it is not unique to a region. 
  ## (two regions can have the same n_i)

  ## May 28
  ## Add Jacobian of the dynamical system to assess the stability of the equilibrium points 
  
  ### Jacobian using original formulation without area size 
  jacobian_balance <- function(n, par) {
    rho_w <- exp(par[1])
    rho_d <- exp(par[2])
    rho_e <- exp(par[3])
    
    wD      <- par[4]
    dD_dest <- par[5]
    dD_src  <- par[6]
    eD      <- par[7]
    
    K  <- length(n)
    L  <- log(n) + 1
    L2 <- log(n + 1) + 1
    
    J <- matrix(0, nrow = K, ncol = K)
    
    for (i in seq_len(K)) {
      incoming <- sum(L[-i]^dD_src)
      
      # diagonal: dF_i / dn_i
      J[i, i] <- rho_w * wD * L[i]^(wD - 1) / n[i] +
        rho_d * dD_dest * L2[i]^(dD_dest - 1) * incoming / (n[i] + 1) -
        rho_e * eD * L[i]^(eD - 1) / n[i]
      
      # off-diagonal: dF_i / dn_k, k != i
      for (k in seq_len(K)) {
        if (k != i) {
          J[i, k] <- rho_d * L2[i]^dD_dest *
            dD_src * L[k]^(dD_src - 1) / n[k]
        }
      }
    }
    
    J
  }

### Given target n*, find the parameters numerically using the original formulation 
# target_n <- c(20, 40, 80)
  target_n <- c(4, 8, 16)
# target_n <- c(40, 20, 80)
# target_n <- c(20, 20, 20)

# base rates in par to be optimized  are in log scale, 
# so optimization does not give negative base rates
objective <- function(par) {
  rho_w <- exp(par[1])
  rho_d <- exp(par[2])
  rho_e <- exp(par[3])
  
  wD <- par[4]
  dD_dest <- par[5]
  dD_src <- par[6]
  eD <- par[7]
  
  L <- log(target_n) + 1
  L_2 <- log(target_n+1) + 1
  K <- length(target_n)
  
  balance <- numeric(K)
  
  for (i in seq_len(K)) {
    incoming <- sum(L[-i]^dD_src)
    
    balance[i] <- rho_w * L[i]^wD +
      rho_d * L_2[i]^dD_dest * incoming -
      rho_e * L[i]^eD
  }
  # optimize balance^2 function because if not, then for example
  # balance in reg 1=-10, balance in reg 2= 10, their sum = 0 even tho we want 
  # their sum = 0 because balance in 1 = 0 and balance in 2 = 0. 
  sum(balance^2) 
}

fit <- optim(
  par = c(log(0.01), log(0.001), log(0.01),
          -2, -1, -1.5, 1),
  fn = objective,
  method = "L-BFGS-B",
  lower = c(log(1e-6), log(1e-6), log(1e-6), -5, -5, -5, 0),
  upper = c(log(1),    log(1),    log(1),     5,  5,  5,  5),
  control = list(maxit = 10000)
)

J <- jacobian_balance(target_n, fit$par)
eigen(J)$values
all(Re(eigen(J)$values) < 0)



## Adding a paramaeter with value unique to each region e.g., area size
## helps solving the exchangeability issue in the case where regions have
## unequal equilibrium diversity 

# target_n <- c(8, 16, 32)
target_n <- c(5, 10, 15)
# target_n <- c(4, 8, 16)
# area     = c(10/100,50/100,100/100) #use relative area size to avoid negative log
area     = c(0.01,0.05,1) #use relative area size to avoid negative log

# base rates in par to be optimized  are in log scale, 
# so optimization does not give negative base rates

objective_area <- function(par) {
  rho_w = exp(par[1])
  rho_d = exp(par[2])
  rho_e = exp(par[3])
  
  wD      = par[4]
  dD_dest = par[5]
  dD_src  = par[6]
  eD      = par[7]
  
  # diversity per unit area
  n_density = target_n / area
  
  L   = log(n_density) + 1
  L_2 = log(n_density + 1) + 1
  K   = length(target_n)
  
  balance = numeric(K)
  
  for (i in seq_len(K)) {
    incoming = sum(L[-i]^dD_src)
    
    balance[i] <- rho_w * L[i]^wD +
      rho_d * L_2[i]^dD_dest * incoming -
      rho_e * L[i]^eD
  }
  
  sum(balance^2)
}


# fit_area <- optim(
#   par = c(log(0.01), log(0.001), log(0.01),
#           -2, -1, -1.5, 1),
#   fn = objective_area,
#   method = "L-BFGS-B",
#   lower = c(log(1e-6), log(1e-6), log(1e-6), -5, -5, -5, -5),
#   upper = c(log(1),    log(1),    log(1),     5,  5,  5,  5),
#   control = list(maxit = 10000)
# )

fit_area <- optim(
  par = c(log(0.01), log(0.001), log(0.01),
          -2, -1, -1.5, 1),
  fn = objective_area,
  method = "L-BFGS-B",
  lower = c(log(1e-6), log(1e-6), log(1e-6), -5, -5, -5, 0),
  upper = c(log(2),    log(2),    log(2),     0,  0,  0,  5),
  control = list(maxit = 10000)
)
fit_area$par

# to avoid printing the truncated version by default
print(fit_area$par, digits = 22)

# Jacobian using new formulation with area size 
jacobian_area <- function(n, par,area) {
  rho_w <- exp(par[1])
  rho_d <- exp(par[2])
  rho_e <- exp(par[3])
  
  wD      <- par[4]
  dD_dest <- par[5]
  dD_src  <- par[6]
  eD      <- par[7]
  
  # diversity per unit area
  n_density = n / area
  
  K  <- length(n)
  L  <- log(n_density) + 1
  L2 <- log(n_density + 1) + 1
  
  J <- matrix(0, nrow = K, ncol = K)
  
  for (i in seq_len(K)) {
    incoming <- sum(L[-i]^dD_src)
    
    # diagonal: dF_i / dn_i
    J[i, i] <- rho_w * wD * L[i]^(wD - 1) / n[i] +
      rho_d * dD_dest * L2[i]^(dD_dest - 1) * incoming / (n[i] + 1) -
      rho_e * eD * L[i]^(eD - 1) / n[i]
    
    # off-diagonal: dF_i / dn_k, k != i
    for (k in seq_len(K)) {
      if (k != i) {
        J[i, k] <- rho_d * L2[i]^dD_dest *
          dD_src * L[k]^(dD_src - 1) / n[k]
      }
    }
  }
  
  J
}

J <- jacobian_area(target_n, fit_area$par, area)
eigen(J)$values
max(eigen(J)$values)
