library(expm)
#Thsi function returns the smallest singular value s_hat

singval <- function(Psi, B){
  N <- dim(B)[1]
  C <- dim(B)[2]
  conjtr_B<-t(B)
  Gb <- conjtr_B %*% B/N
  Gb1 <- solve(Gb)
  S <- t(Psi) %*% B / N
  # D= Psi'*B(B'B)B'Psi
  D <- S %*% Gb1 %*% t(S)
  
  Gpsi <- t(Psi) %*% Psi / N
  l <- min(eigen(Gpsi)$values)
  So <- solve(sqrtm(Gpsi)) %*% S %*% solve(sqrtm(Gb))
  minSin<- min(svd(So)$d)
  return(minSin)
}

singval(Psi,B)