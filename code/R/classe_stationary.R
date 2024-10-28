stationary.freq.bisse <- function(pars) {
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  
  g <- (lambda0 - mu0) - (lambda1 - mu1)
  eps <- (lambda0 + mu0 + lambda1 + mu1)*1e-14
  if ( abs(g) < eps ) {
    if ( q01 + q10 == 0 )
      p <- 0.5
    else
      p <- q10/(q01 + q10)
  } else {
    roots <- quad_roots(g, q10+q01-g, -q10)
    roots <- roots[roots >= 0 & roots <= 1]
    if ( length(roots) > 1 )
      p <- NA
    else
      p <- roots
  }
  c(p, 1-p)
}
########################################
stationary.freq.classe <- function(pars, k) {
  if (k == 2) {
    g <- (sum(pars[1:3]) - pars[7]) - (sum(pars[4:6]) - pars[8])
    eps <- sum(pars[1:8]) * 1e-14
    ss1 <- pars[9]  + 2*pars[3] + pars[2]  # shift from 1
    ss2 <- pars[10] + 2*pars[4] + pars[5]  # shift from 2
    
    if ( abs(g) < eps ) {
      if (ss1 + ss2 == 0) 
        eqfreq <- 0.5
      else
        eqfreq <- ss2/(ss1 + ss2)
      eqfreq <- c(eqfreq, 1 - eqfreq)
    } else {
      roots <- quad_roots(g, ss2 + ss1 - g, -ss2)
      eqfreq <- roots[roots >= 0 & roots <= 1]
      if ( length(eqfreq) > 1 )
        eqfreq <- NA
      else
        eqfreq <- c(eqfreq, 1 - eqfreq)
    }
  } else { ## also works for k=2, but much slower
    eqfreq <- stationary.freq.classe.ev(pars, k)
  }
  eqfreq
}

stationary.freq.classe.ev <- function(pars, k) {
  A <- projection.matrix.classe(pars, k)
  ## continuous time, so the dominant eigenvalue is the largest one
  ## return its eigenvector, normalized
  evA <- eigen(A)
  i <- which(evA$values == max(evA$values))
  evA$vectors[,i] / sum(evA$vectors[,i])
}

projection.matrix.classe <- function(pars, k) {
  A <- matrix(0, nrow=k, ncol=k)
  
  nsum <- k*(k+1)/2
  kseq <- seq_len(k)
  pars.lam <- pars[seq(1, nsum*k)]
  pars.mu <- pars[seq(nsum*k+1, (nsum+1)*k)]
  pars.q <- pars[seq((nsum+1)*k+1, length(pars))]
  
  ## array indices of lambda's in parameter vector
  idx.lam <- cbind(rep(kseq, each=nsum), rep(rep(kseq, times=seq(k,1,-1)), k),
                   unlist(lapply(kseq, function(i) i:k)))
  ## transpose of matrix indices of q's in parameter vector
  idx.q <- cbind(unlist(lapply(kseq, function(i) (kseq)[-i])), 
                 rep(kseq, each=k-1))
  
  ## take care of off-diagonal elements
  for (n in seq_len(nsum*k)) {
    ## add this lambda to A[daughter states, parent state]
    ## (separate steps in case the daughter states are the same)
    r <- idx.lam[n,]
    A[r[2], r[1]] <- A[r[2], r[1]] + pars.lam[n]
    A[r[3], r[1]] <- A[r[3], r[1]] + pars.lam[n]
  }
  A[idx.q] <- A[idx.q] + pars.q
  
  ## fix the diagonal elements
  diag(A) <- 0
  diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) 
    sum(pars.lam[seq((i-1)*nsum+1, i*nsum)]) - pars.mu[i]))
  A
}
###########################################################
stationary.freq.geosse <- function(pars) {
  sA  <- pars[1]
  sB  <- pars[2]
  sAB <- pars[3]
  xA  <- pars[4]
  xB  <- pars[5]
  dA  <- pars[6]
  dB  <- pars[7]
  
  A <- matrix(c(
    -xA - xB - sAB,   dA,             dB,
    sA + xB + sAB,   sA - dA - xA,   0,
    sB + xA + sAB,   0,              sB - dB - xB
  ), nrow=3, byrow=TRUE)
  
  # A <- matrix(c(
  #   -0.5*sA - 0.5*sB -xA - xB - sAB,   dA,             dB,
  #   0.5*sA + xB + sAB,   sA - dA - xA,   0,
  #   0.5*sB + xA + sAB,   0,              sB - dB - xB
  # ), nrow=3, byrow=TRUE)
  
  # A <- matrix(c(
  #   -sA-sB-xA - xB - sAB,   dA,             dB,
  #   sA + xB + sAB,   sA - dA - xA,   0,
  #   sB + xA + sAB,   0,              sB - dB - xB
  # ), nrow=3, byrow=TRUE)
  # A <- matrix(c(
  #   sA+sB-xA - xB - sAB,   dA,             dB,
  #    xB + sAB,   - dA - xA,   0,
  #    xA + sAB,   0,              - dB - xB
  # ), nrow=3, byrow=TRUE)
  ## continuous time, so the dominant eigenvalue is the largest one
  ## return its eigenvector, normalized  
  evA <- eigen(A)
  i <- which(evA$values == max(evA$values))
  evA$vectors[,i] / sum(evA$vectors[,i])
}
