library(mgcv)
library(parallel)

TSCAN_group <- function(expr,pseudotime,branch, parallel = FALSE, n.cores = 8) {
  if (parallel){
    pval <- mclapply(row.names(expr),function(g) {
      x <- expr[g,]
      m1 <- mgcv::gam(x~s(pseudotime,k=3,by=branch))
      m2 <- mgcv::gam(x~s(pseudotime,k=3))
      c(pchisq(m2$deviance - m1$deviance, m2$df.residual - m1$df.residual,lower.tail = F),m2$deviance - m1$deviance)
    },mc.cores = n.cores)
  } else {
    pval <- lapply(row.names(expr),function(g) {
      x <- expr[g,]
      m1 <- mgcv::gam(x~s(pseudotime,k=3,by=branch))
      m2 <- mgcv::gam(x~s(pseudotime,k=3))
      c(pchisq(m2$deviance - m1$deviance, m2$df.residual - m1$df.residual,lower.tail = F),m2$deviance - m1$deviance)
    })
  }
  names(pval) <- row.names(expr)
  res <- do.call(rbind,pval)
  res <- data.frame(pval=res[,1],chisq=res[,2])
  res$fdr <- p.adjust(res[,1],method='fdr')
  res
}

TSCAN_time <- function(expr, pseudotime, parallel = FALSE, n.cores = 8) {
  if (parallel){
    pval <- mclapply(row.names(expr),function(g) {
      print(g)
      x <- expr[g,]
      model <- mgcv::gam(x~s(pseudotime,k=3))
      c(pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F),model$null.deviance - model$deviance)
    },mc.cores=detectCores())
  } else {
    pval <- lapply(row.names(expr),function(g) {
      x <- expr[g,]
      model <- mgcv::gam(x~s(pseudotime,k=3))
      c(pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F),model$null.deviance - model$deviance)
    })
  }
  names(pval) <- row.names(expr)
  res <- do.call(rbind,pval)
  res <- data.frame(pval=res[,1],chisq=res[,2])
  res$fdr <- p.adjust(res[,1],method='fdr')
  res
}

