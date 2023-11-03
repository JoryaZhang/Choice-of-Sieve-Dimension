#P_K:B * inv(B′*B)*t(B)
P_K<- function(B){
  BtB <- t(B) %*% B
  # Compute the inverse of (B_K^T * B_K)
  inv_BtB <- solve(BtB)
  # Compute P_K
  P_K <- B %*%inv_BtB %*% t(B)
  return(P_K)
}


M_J <- function(Psi, B) {
  # Compute Psi_J^T * P_K_J * Psi_J
  P_K <- P_K(B)
  mat1 <- t(Psi) %*% P_K %*% Psi
  # Compute the inverse of Psi_J^T * P_K_J * Psi_J
  inv_mat1 <- solve(mat1)
  # Compute M_J
  M_J <- inv_mat1 %*% t(Psi) %*% P_K
  return(M_J)
}


#h_hat = psi(x)*M_j*Y
h_hat <- function(Psi, M_J, Y,n) {
  h_hat <- Psi[n,] %*% M_J %*% Y
  return(h_hat)
}
#u=Y_i-h_hat(Xi)
u_J<- function(Psi, M_J, Y){
  ulis <-c()
  n <- dim(Psi)[1]
  for (i in 1:n) {
    h<-h_hat(Psi,M_J,Y,i)
    u <- Y[i]-h
    ulis <- c(ulis, u)
  }
  return(ulis)
}


#This generate the multiplier bootsrtap version of u^J
#vp is N(0,1).
multiboot_u <- function(Psi, M_J, Y, vp) {
  u_J <- u_J(Psi, M_J, Y)
  n <- length(u_J)
  
  # Compute u^*(J)
  u_star_J <- u_J * vp
  
  return(u_star_J)
}


U_JJ2 <- function(u_J, u_J2) {
  n <- length(u_J)
  diagonal_entries <- sapply(1:n, function(i) u_J[[i]] * u_J2[[i]])
  U_JJ2 <- diag(diagonal_entries)
  
  return(U_JJ2)
}


# Compute D_J
D_J<- function(Psi, Y, B, n){
  P_K <- P_K(B)
  M_J <- M_J(Psi,B)
  u_J <- u_J(Psi,M_J,Y)
  D_J <- Psi[n,] %*% M_J %*% u_J
  return(D_J)
}

#Compute D_J_star
D_J_star <- function(Psi, Y, B, n, vp){
  P_K <- P_K(B)
  M_J <- M_J(Psi,B)
  u_star <- multiboot_u(Psi, M_J, Y, vp)
  D_J <- Psi[n,] %*% M_J %*% u_star
  return(D_J)

}



  
