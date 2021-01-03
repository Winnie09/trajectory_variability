difftest <- function (data, TSCANorder, df = 3) {
  # data: gene by cell log-normalized gene expression
  # TSCANorder: a character vector of cells ordered by pseudotime in a branch
  ptime <- 1:length(TSCANorder)
  data <- data[, TSCANorder]
  pval <- mclapply(rownames(data), function(x) {
    x <- data[x,]
    model <- mgcv::gam(x ~ s(ptime, k = 3))
    pchisq(model$null.deviance - model$deviance, model$df.null - 
             model$df.residual, lower.tail = F)
  },mc.cores=detectCores())
  names(pval) <- rownames(data)
  pval <- unlist(pval)
  qval <- p.adjust(pval, method = "fdr")
  data.frame(pval = pval, qval = qval)
}
