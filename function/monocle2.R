library(parallel)
monocle2_group <- function(expr, pseudotime, branch, parallel = FALSE, n.cores = 8) {
  library(VGAM)
  if (parallel){
    pval <- mclapply(row.names(expr),function(g) {
      x <- expr[g,]
      m1 <- suppressWarnings(VGAM::vglm(x~sm.ns(pseudotime, df = 3)*branch, epsilon=1e-1, family='uninormal'))
      m2 <- suppressWarnings(VGAM::vglm(x~sm.ns(pseudotime, df = 3), epsilon=1e-1, family='uninormal'))  
      lrt <- VGAM::lrtest(m1,m2) 
      c(lrt@Body["Pr(>Chisq)"][2,],lrt@Body["Chisq"][2,])
    },mc.cores = n.cores)
  } else {
    pval <- lapply(row.names(expr),function(g) {
      x <- expr[g,]
      m1 <- suppressWarnings(VGAM::vglm(x~sm.ns(pseudotime, df = 3)*branch, epsilon=1e-1, family='uninormal'))
      m2 <- suppressWarnings(VGAM::vglm(x~sm.ns(pseudotime, df = 3), epsilon=1e-1, family='uninormal'))  
      lrt <- VGAM::lrtest(m1,m2) 
      c(lrt@Body["Pr(>Chisq)"][2,],lrt@Body["Chisq"][2,])
    })
  }
  names(pval) <- row.names(expr)
  res <- do.call(rbind,pval)
  res <- data.frame(pval=res[,1],chisq=res[,2])
  res$fdr <- p.adjust(res[,1],method='fdr')
  res
}

monocle2_time <- function(expr,pseudotime) {
  pval <- mclapply(row.names(expr),function(g) {
    x <- expr[g,]
    m1 <- suppressWarnings(VGAM::vglm(x~sm.ns(pseudotime, df=3), epsilon=1e-1, family='uninormal'))
    m2 <- suppressWarnings(VGAM::vglm(x~1, epsilon=1e-1, family='uninormal'))  
    lrt <- VGAM::lrtest(m1,m2) 
    c(lrt@Body["Pr(>Chisq)"][2,],lrt@Body["Chisq"][2,])
  },mc.cores=8)
  names(pval) <- row.names(expr)
  res <- do.call(rbind,pval)
  res <- data.frame(pval=res[,1],chisq=res[,2])
  res$fdr <- p.adjust(res[,1],method='fdr')
  res
}

