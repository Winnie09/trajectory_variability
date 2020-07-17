testpt_Variable <- function(expr, cellanno, pseudotime, design, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=10) {
  # pseudotime: a numeric vector of pseudotime, and the names are cell names
  set.seed(12345)
  fitfunc <- function(iter) {
    if (iter==1) {
      fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)  
    } else {
      dn <- paste0(as.vector(design),collapse = '_')
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
      fitpt(expr=perexpr, cellanno=percellanno, pseudotime=psn, design=perdesign, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)
    }
  }
  
  if (ncores == 1){
    fit <- lapply(1:(permuiter+1),fitfunc)
  } else {
    fit <- mclapply(1:(permuiter+1),fitfunc, mc.cores = ncores)
  }
  
  orifit <- fit[[1]]
  knotnum <- orifit$knotnum
  orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)[row.names(expr)]
  perll <- sapply(2:(permuiter+1),function(i) sapply(fit[[i]]$parameter,function(i) unname(i$ll),USE.NAMES = F)[row.names(expr)])

  pval <- sapply(1:nrow(perll), function(i) {
    z <- perll[i,]
    den <- density(z)$bw
    mean(pnorm(orill[i], z, sd=den,lower.tail = F))
  })
  fdr <- p.adjust(pval,method='fdr')
  names(fdr) <- row.names(perll)
  foldchange <- orill - rowMeans(perll)
  return(list(fdr = fdr, orill=orill, perll = perll, knotnum = knotnum, foldchange = foldchange, parameter=orifit$parameter, pseudotime = pseudotime))
}


