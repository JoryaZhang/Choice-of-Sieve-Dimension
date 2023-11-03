#In the simulations and empirical application we have a scalar X(d=1)
#use a cubic B-spline basis (i.e., r = 4)
#We use a cubic B-spline basis (r = 4) for J and quartic(r=5) B-spline basis for K
#K is pinned down by J
J_max<-function(r,d,X,W,rJ,rK){
  # s <- shat(Psi,Y,B)$shat
  Tau <- Tset(r,d)
  J <- Tau$J
  K <- Tau$K
  i=1
  jlis<-c()
  while (TRUE) {
    if(i+1>length(J)){
      break
    }
    j1=J[i]
    j2=J[i+1]
    k1=K[i]
    k2=K[i+1]
    SJ <- j1-rJ
    SK <- k1-rK
    SJ_m <- j2-rJ
    SK_m <- k2-rK
    #generate matrices
    Psi <- bsplinemat_quantile(X,SJ,rJ)$XX
    B <- bsplinemat_quantile(W,SK,rK)$XX
    Psi_m <- bsplinemat_quantile(X,SJ_m,rJ)$XX
    B_m <- bsplinemat_quantile(W,SK_m,rK)$XX
    #compute shat and shat+
    shat <- singval(Psi,B)
    shat_m <- singval(Psi_m, B_m)
    N <- dim(B)[1]
    lhs <- j1*sqrt(log(j1))*(1/shat)
    rhs <- j2*sqrt(log(j2))*(1/shat_m)
    if (lhs <= 10*sqrt(N) && rhs > 10*sqrt(N)){
      jlis <- c(jlis, j1)
    }
    i <- i+1
  }
  return(min(jlis))
}

#This function returns TJ_hat set
TJ_hat <- function(r,d,J_max) {
  log_J_max <- log(J_max)
  lb <- 0.1 * (log_J_max)^2
  Jset<-Tset(r,d)$J
  Kset <- Tset(r,d)$K
  l <- length(Jset)
  j_hat<-c()
  k_hat <- c()
  for (i in 1:l){
    if(Jset[i]>=lb && Jset[i]<= J_max){
      j_hat<-c(j_hat,Jset[i])
      k_hat <- c(k_hat,Kset[i])
    }
  }
  return(list("j_hat"=j_hat, "k_hat"=k_hat))
}


