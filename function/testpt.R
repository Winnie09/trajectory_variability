testpt <- function(expr, cellanno, pseudotime, design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores()) {
  orifit <- fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores)
  knotnum <- orifit$knotnum
  orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)
  orill <- orill[row.names(expr)]
  dn <- paste0(as.vector(design),collapse = '_')
  perll <- lapply(1:permuiter,function(did) {
    print(did)
    perdn <- dn
    while(perdn==dn) {
      perid <- sample(1:nrow(design))
      perdesign <- design[perid,,drop=F]
      perdn <- paste0(as.vector(perdesign),collapse = '_')  
    }
    row.names(perdesign) <- row.names(design)
    sampcell <- sample(1:ncol(expr),replace=T)
    perexpr <- expr[,sampcell,drop=F]
    percellanno <- cellanno[sampcell,,drop=F]
    psn <- pseudotime[colnames(expr)]
    psn <- psn[sampcell]
    colnames(perexpr) <- percellanno[,1] <- names(psn) <- paste0('cell_',1:length(psn))
    perfit <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=psn, design=perdesign, knotnum=knotnum, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)
    perfit <- sapply(perfit$parameter,function(i) unname(i$ll))
    perfit <- perfit[row.names(expr)]
  })
  perll <- do.call(cbind,perll)
  pval <- sapply(1:nrow(perll), function(i) pnorm(orill[i],mean(perll[i,]),sd(perll[i,]),lower.tail = F))
  fdr <- p.adjust(pval,method='fdr')
  names(fdr) <- row.names(perll)
  foldchange <- orill - rowMeans(perll)
  return(list(fdr = fdr, perll = perll, knotnum = knotnum, foldchange = foldchange, parameter=orifit$parameter))
}

