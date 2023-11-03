#This function computes the sup statistics of each independent draw of varpi
#J set and x are settled
sigmahat_sq <- function(Psi, Y, B, n){
  psi <- Psi[n,]
  M_J <- M_J(Psi, B)
  u_J <- u_J(Psi,M_J,Y)
  U_JJ2 <- U_JJ2(u_J, u_J)
  s <- t(psi) %*% M_J %*% U_JJ2 %*% t(M_J) %*% psi
  return(s)
}

sigmatilda <- function(Psi, Psi2, Y, B, B2, n){
  psi <- Psi[n,]
  psi2 <- Psi2[n,]
  M_J <- M_J(Psi, B)
  M_J2 <- M_J(Psi2, B2)
  u_J <- u_J(Psi,M_J,Y)
  u_J2 <- u_J(Psi2,M_J2,Y)
  U_JJ2 <- U_JJ2(u_J, u_J2)
  s <- t(psi) %*% M_J %*% U_JJ2 %*% t(M_J2) %*% psi2
  return(s)
}

sigmajj2 <- function(Psi, Psi2, Y, B, B2, n){
  s1 <- sigmahat_sq(Psi, Y, B, n)
  s2 <- sigmahat_sq(Psi2, Y, B2, n)
  st <- sigmatilda(Psi, Psi2, Y, B, B2, n)
  s <- s1 + s2 - 2*st
  return(sqrt(s))
}

Jmax <- J_max(4,1,X,W,4,5)
Jhat <- TJ_hat(4,1,Jmax)$j_hat
Khat <- TJ_hat(4,1,Jmax)$k_hat
sup_abs <- function(X, W, Y, Jhat, Khat, rJ, rK){
  suplis <- c()
  n <- length(X)
  vp <- rnorm(n)
  for (r in 1:n){
    i=1
    xlis<-c()
    while (TRUE) {
      if(i>=length(Jhat)){
        break
      }
      #fix a j
      j <- Jhat[i]
      k <- Khat[i]
      SJ <- j-rJ
      SK <- k-rK
        #search for all j2>j
      for (u in i:length(Jhat)){
          SJ2 <- Jhat[u]-rJ
          SK2 <- Khat[u]-rK
          Psi <- bsplinemat_quantile(X,SJ,rJ)$XX
          B <- bsplinemat_quantile(W,SK,rK)$XX
          Psi2 <- bsplinemat_quantile(X,SJ2,rJ)$XX
          B2 <- bsplinemat_quantile(W,SK2,rK)$XX
          D_dif <- D_J_star(Psi, Y, B, r, vp)-D_J_star(Psi2, Y, B2, r, vp)
          sigjj2 <- sigmajj2(Psi, Psi2, Y, B, B2, r)
          if (sigjj2 ==0) next
          result <- abs(D_dif / sigjj2)
          xlis <- c(xlis, result)
          }
      i <- i+1
    }
    suplis <- c(suplis, max(xlis))
  }
  return(max(suplis))
}
#It taked too long to run 1027 numbers, so I choose the first 100 rows for test
a <- sup_abs(X[1:100],W[1:100], Y[1:100],Jhat,Khat,4,5)

#This function returns the （1-alpha) quantile of sup statistic across t draws
#t: test times (also t times of draws)
#alpha: min(0.5, sqrt(log(Jmax)/Jmax))
alpha <- min(0.5, sqrt(log(Jmax)/Jmax))
statquantile <- function(alpha, t, X, W, Y, Jhat, Khat, rJ, rK){
  data<-c()
  for (i in 1:t) {
    data <- c(data, sup_abs(X, W, Y, Jhat, Khat, rJ, rK) )
  }
  theta <- quantile(data, 1 - alpha)
  return(theta)
}
t<-100
theta <- statquantile(1-alpha, t, X[1:100],W[1:100], Y[1:100],Jhat,Khat,4,5)
print(theta)
