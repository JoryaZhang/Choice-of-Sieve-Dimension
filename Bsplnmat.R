bsplinemat <- function(x, m, r, kts){
  m <- m+2
  XX <- bsplnbasis(x, m, r, kts)$XX
  DX <- bsplnbasis(x, m, r, kts)$DX
  DDX <- bsplnbasis(x, m, r, kts)$DDX
  return(list("XX"=XX, "DX"=DX, "DDX"=DDX))
}
bsplinemat_quantile <- function(x, m, r) {
  if (m == 0) {
    kts <- c(rep(0, r), rep(1, r))
  } else {
    quantiles <- quantile(x, (1:1:m) / (m + 1))
    kts <- c(rep(0, r), quantiles, rep(1, r))
  }
  result <- bsplinemat(x, m, r, kts)
  list("XX" = result$XX, "DX" = result$DX, "DDX" = result$DDX, "kts" = kts)
}

# Call the bsplinemat_quantile function
# result <- bsplinemat_quantile(x, m, r)
# XX <- result$XX
# DX <- result$DX
# DDX <- result$DDX
# kts <- result$kts
