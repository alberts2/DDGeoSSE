#Example:
# using parameters derived from lemma 1
# should give us ~100 species at equlibrium per region
# solve_balance(1.741029,0.01,0.01,-1,2,0,-1,num_reg=3,c(1,100),tol = 1e-6,max_expand = 2)
#
# using parameters derived from lemma 2
# should give us ~100 species at equlibrium per region
# solve_balance(12.47404,0.04,0.03,-0.5,3,0,-0.5,num_reg=3,c(1,100),tol = 1e-6,max_expand = 2)

# note that this script assumes that no information known about n_j*. So, we always assume 
# that n_i* = n_j* for all i and j. 

# note if you get Error in if (fa * fb < 0) { : missing value where TRUE/FALSE needed
# that means the lower bound becomes too small such that we have log(0)

# 8 June: added parameters from benchmarking test with 3 to 7 regions 
num_region = 3
solve_balance(1.0,0.1,0.2,-0.6,0.5,0,-0.5,num_region,c(1,100),tol = 1e-6,max_expand = 2)$root*num_region
solve_balance <- function(rho_w, rho_d, rho_e,
                          w_d, e_d, d_d_src, d_d_dest,
                          num_reg,
                          interval,
                          tol,
                          max_expand){
  balance <- function(n) {
    left_eq = rho_w*(log(n)+1)^(w_d) + (num_reg-1)*rho_d*(log(n)+1)^(d_d_src)*(log(n+1)+1)^(d_d_dest)
    #
    right_eq = rho_e*(log(n)+1)^(e_d)
    #
    left_eq - right_eq
  }
  #
  a <- interval[1]
  b <- interval[2]
  
  # Try expanding until sign change is found
  for (i in 1:max_expand) {
    fa <- balance(a)
    fb <- balance(b)
    
    # If signs differ, uniroot can find an exact zero
    if (fa * fb < 0) {
      root <- uniroot(balance, c(a, b), tol = tol)$root
      return(list(
        root = root,
        difference = abs(balance(root)),
        interval = c(a, b),
        expanded = i - 1,
        status = "root found"
      ))
    }
    # if no sign change, expand interval
    a <- a / 2
    b <- b * 2
  }
  # fallback: interval expansion failed
  return(list(
    root = NA,
    difference = NA,
    interval = c(a, b),
    expanded = max_expand,
    status = "no root (no sign change detected)- try changing "
  ))
}