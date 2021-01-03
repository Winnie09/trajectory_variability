get_gam_fit <- function(expr, TSCANorder, n.cores = 4){
  ptime <- 1:length(TSCANorder)
  data <- expr[, TSCANorder]
  tmp <- mclapply(rownames(data), function(x) {
    x <- data[x,]
    model <- mgcv::gam(x ~ s(ptime, k = 3))
    fitted(model)
  }, mc.cores = n.cores)
  tmp <- do.call(rbind, tmp)
  dimnames(tmp) <- dimnames(data)
  tmp
}

