scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- apply(data,1,sd)
  (data - cm) / csd
}
