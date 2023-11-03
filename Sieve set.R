#This function defines the sieve dimension set with rth order and d-dimensional
#It returns J /in T and K is pinned down by J
#K is pinned down by J
#l_w=⌈l+q⌉, where l is the resolution level of J, here we pick q=1
Tset <- function(r, d) {
  # Initialize an empty vector to store the elements of T
  lis <- c()
  Klis <- c()
  l <- 0
  while (TRUE) {
    J <- (2^l + r - 1)^d
    K <- (2^(l+1)+r)^d
    if (J > 100) {
      break
    }
    # Append the element to the vector T
    lis <- c(lis, J)
    Klis <- c(Klis, K)
    l <- l + 1
  }
  return(list("J"=lis, "K"=Klis))
}