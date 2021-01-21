library(mgcv)
library(parallel)

TSCAN_group <- function(expr,pseudotime,branch, parallel = FALSE, n.cores = 8) {
 # expr: saver imputed matrix
 # count: count matrix
 # pseudotime: numeric vector (1,2,3,4....) with names same as colnames(expr)
 # branch: 0,1 vector indicating whether each cell is from group 1 or 2, can get from as.numeric(sub(':.*','',colnames(expr)) %in% paste0('BM',c(1,2,5,6)))
 # cell_coords: the pca you sent me, only use the first 4 (if correct) dimensions
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

TSCAN_testvar <- function(expr, cellanno, design, pseudotime){
  ## count: gene by cell count matrix
  ## cellanno: 1st column cell name (chr), 2nc column sample name (chr)
  ## design: 1st column the group division (num), colname are the variable name (e.g. 'group'). rownames are sample names.
  ## pseudotime: 1st column cell name (chr), 2nd column pseudotime (num).
  design = cbind(1,design)
  psn <- pseudotime[,2]
  names(psn) <- pseudotime[,1]
  branch <- sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res <- TSCAN_group(expr, psn, branch)
  return(res)
}

