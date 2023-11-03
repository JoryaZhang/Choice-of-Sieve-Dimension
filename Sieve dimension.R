#This function returns sup J satisfied (11)
sup_J <- function(theta, X, W, Y, rJ, rK, r,d ){
  Jmax <- J_max(r,d,X,W,rJ,rK)
  Jhat <- TJ_hat(r,d,Jmax)$j_hat
  Khat <- TJ_hat(r,d,Jmax)$k_hat
  jlis<- c()
  n <- length(X)
  for (r in 1:n){
    i=1
    xlis<-c()
    while (TRUE) {
      if(i+1>length(Jhat)){
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
        D_dif <- D_J(Psi, Y, B, r)-D_J(Psi2, Y, B2, r)
        sigjj2 <- sigmajj2(Psi, Psi2, Y, B, B2, r)
        if (sigjj2==0) next
        result <- abs(D_dif / sigjj2)
        xlis <- c(xlis, result)
        result_max <- max(xlis)
        if(result_max < 1.1*theta){
          jlis <- c(jlis,j)
        }
      }
      
      i <- i+1
    }
  }
  return(min(jlis))
}

Jtilda <- function(theta, X, W, Y, rJ, rK, r,d){
  supJ <- sup_J(theta, X, W, Y, rJ, rK, r,d)
  jmax <- J_max(r,d,X,W,rJ,rK)
  jlis <- Tset(r,d)$J
  Jn <- max(jlis[jlis < jmax])
  return(min(Jn,supJ))
}




