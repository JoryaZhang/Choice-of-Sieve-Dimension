bsplnbasis <- function(x, m, r, kts = NULL) {
    N=length(x)
    if (is.null(kts)) {
      kts <- numeric(m + 2 * r - 2)#generate m+2r-2 zero vectors
      kts[1:(r - 1)] <- numeric(r - 1) #from position 1 to position r-1, kts=0
      kts[(m + r):(m + 2 * r - 2)] <- numeric(r - 1) + 1
      kts[(r:(m + r - 1))] <- seq(0, 1, length.out = m) #from position r to m+r-1 equally divided [0,1] by m
    }
    # initialize for recursion
    # generate r matrices of 0 which has N rows and m+2r-2 columns
    BB <- array(0, dim = c(N, m + 2 * r - 2, r)) 
    # give the initial lower bound and upper bound for x[i]
    for (i in 1:N) {
      ix <- min(which(x[i] >= kts[(r):(r + m)] & x[i] <= kts[(r + 1):(r + m + 1)]))
      BB[i, ix + r - 1, 1] <- 1
    }
    
    
    # recursion
    for (j in 2:r) {
      for (i in 1:(m + 2 * r - 2 - j)) {
        if ((i + j + 1) <= (m + 2 * r)) {
          if ((kts[(i + j - 1)] - kts[i]) != 0) {
            a1 <- (x - kts[i]) / (kts[(i + j - 1)] - kts[i])
          } else {
            a1 <- numeric(N)
          }
          if ((kts[(i + j)] - kts[(i + 1)]) != 0) {
            a2 <- (x - kts[(i + 1)]) / (kts[(i + j)] - kts[(i + 1)])
          } else {
            a2 <- numeric(N)
          }
          BB[, i, j] <- a1 * BB[, i, j - 1] + (1 - a2) * BB[, (i + 1), (j - 1)]
        } else if ((i + j) <= (m + 2 * r)) {
          if ((kts[(i + j)] - kts[i]) != 0) {
            a1 <- (x - kts[i]) / (kts[(i + j)] - kts[i])
          } else {
            a1 <- numeric(N)
          }
          BB[, i, j] <- a1 * BB[, i, j - 1]
        }
      }
    }
    XX <- BB[, 1:(m + r - 2), r]
    
# print(BB)
# XX <- BB[, 1:(m + r - 2), r]
# print(XX)

  DX <- matrix(0, nrow = N, ncol = m + r - 2)
  
  if (r > 1) {
    for (i in 1:(m + r - 2)) {
      if ((kts[(i + r - 1)] - kts[i]) != 0) {
        a1 <- rep(1, N) / (kts[(i + r - 1)] - kts[i])
      } else {
        a1 <- numeric(N)
      }
      if ((kts[(i + r)] - kts[(i + 1)]) != 0) {
        a2 <- rep(1, N) / (kts[(i + r)] - kts[(i + 1)])
      } else {
        a2 <- numeric(N)
      }
      
      if (i < (m + r)) {
        DX[, i] <- (r - 1) * (a1 * BB[, i, (r - 1)] - a2 * BB[, (i + 1), (r - 1)])
      } else {
        DX[, i] <- (r - 1) * (a1 * BB[, i, (r - 1)])
      }
    }
  }
  
  DDX <- matrix(0, nrow = N, ncol = m + r - 2)
  
  if (r > 2) {
    for (i in 1:(m + r - 2)) {
      if ((kts[(i + r - 1)] - kts[i]) != 0) {
        a1 <- rep(1, N) / (kts[(i + r - 1)] - kts[i])
      } else {
        a1 <- numeric(N)
      }
      if ((kts[(i + r)] - kts[(i + 1)]) != 0) {
        a2 <- rep(1, N) / (kts[(i + r)] - kts[(i + 1)])
      } else {
        a2 <- numeric(N)
      }
      if ((kts[(i + r - 1)] - kts[(i + 1)]) != 0) {
        a3 <- rep(1, N) / (kts[(i + r - 1)] - kts[(i + 1)])
      } else {
        a3 <- numeric(N)
      }
      if ((kts[(i + r - 2)] - kts[i]) != 0) {
        a4 <- rep(1, N) / (kts[(i + r - 2)] - kts[i])
      } else {
        a4 <- numeric(N)
      }
      if ((kts[(i + r)] - kts[(i + 2)]) != 0) {
        a5 <- rep(1, N) / (kts[(i + r)] - kts[(i + 2)])
      } else {
        a5 <- numeric(N)
      }
      
      DDX[, i] <- (r - 1) * (r - 2) * (a4 * a1 * BB[, i, (r - 2)] + a5 * a2 * BB[, (i + 2), (r - 2)] - (a1 + a2) * a3 * BB[, (i + 1), (r - 2)])
    }
  }
  out<-list("XX"=XX, "DX"=DX, "DDX"=DDX)
  return(out)
}

