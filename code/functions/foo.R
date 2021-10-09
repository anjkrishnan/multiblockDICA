## foo ===================
## count the number of each level
## =======================
foo <- function(x) {
  xuniq <- unique(x)
  N <- length(xuniq)
  res <- rep(NA, N)
  for (i in 1:N) {
    res[i] <- sum(x == xuniq[i])
  }
  return(res)
}
