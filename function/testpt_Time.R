testpt_Time <- function(expr, cellanno, pseudotime, design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores()) {
  design = as.matrix(design[,1,drop=F])
  design[,1] <- 1
  orifit <- fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores,parallel=T)
  knotnum <- orifit$knotnum
  orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)
  orill <- orill[row.names(expr)]
  if (ncores == 1){
    perll <- lapply(1:permuiter,function(did) {
      sampcell <- sample(1:ncol(expr),replace=T)
      perexpr <- expr[,sampcell,drop=F]
      percellanno <- cellanno[sampcell,,drop=F]
      psn <- pseudotime[colnames(expr)]
      psn <- psn[sampcell]
      colnames(perexpr) <- percellanno[,1] <- names(psn) <- paste0('cell_',1:length(psn))
      prepsn <- psn
      set.seed(did)
      names(prepsn) <- sample(names(prepsn))
      perfit <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=prepsn, design=design, knotnum=NULL, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)
      perfit <- sapply(perfit$parameter,function(i) unname(i$ll))
      perfit <- perfit[row.names(expr)]
    })
  } else {
    perll <- mclapply(1:permuiter,function(did) {
      sampcell <- sample(1:ncol(expr),replace=T)
      perexpr <- expr[,sampcell,drop=F]
      percellanno <- cellanno[sampcell,,drop=F]
      psn <- pseudotime[colnames(expr)]
      psn <- psn[sampcell]
      colnames(perexpr) <- percellanno[,1] <- names(psn) <- paste0('cell_',1:length(psn))
      prepsn <- psn
      set.seed(did)
      names(prepsn) <- sample(names(prepsn))
      perfit <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=prepsn, design=design, knotnum=NULL, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)
      perfit <- sapply(perfit$parameter,function(i) unname(i$ll))
      perfit <- perfit[row.names(expr)]
    }, mc.cores = ncores)
  }
  perll <- do.call(cbind,perll)
  pval <- sapply(1:nrow(perll), function(i) pnorm(orill[i],mean(perll[i,]),sd(perll[i,]),lower.tail = F))
  fdr <- p.adjust(pval,method='fdr')
  names(fdr) <- row.names(perll)
  foldchange <- orill - rowMeans(perll)
  return(list(fdr = fdr, perll = perll, knotnum = knotnum, foldchange = foldchange, parameter=orifit$parameter))
}


