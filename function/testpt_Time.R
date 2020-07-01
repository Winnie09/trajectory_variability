testpt_Time <- function(expr, cellanno, pseudotime, design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=10) {
  expr <- expr[, names(pseudotime)]
  cellanno <- cellanno[match(names(pseudotime), cellanno[,1]), ]
  design = as.matrix(design[,1,drop=F])
  design[,1] <- 1
  #psn <- pseudotime[,2]
  #names(psn) <- pseudotime[,1]
  #pseudotime <- psn
  set.seed(12345)
  orifit <- fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores,parallel=T)
  knotnum <- orifit$knotnum
  orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)
  orill <- orill[row.names(expr)]
  if (ncores == 1){
    perll <- lapply(1:permuiter,function(did) {
      ### permute time within each sample , 20200630
      perpsn <- sapply(rownames(design), function(s){
        tmpid <- cellanno[cellanno[,2] == s, 1]  # subset cells
        tmppsn <- pseudotime[names(pseudotime) %in% tmpid] # subset time
        names(tmppsn) <- sample(names(tmppsn)) # permute time
        tmppsn
      })
      names(perpsn) <- NULL
      perpsn <- unlist(perpsn)
      perpsn <- perpsn[colnames(expr)]
      ### bootstrap within each sample , 20200630
      sampcell <- as.vector(unlist(lapply(unique(cellanno[,2]), function(p){
        sample(which(cellanno[,2] == p), replace = T)
      })))
      perexpr <- expr[,sampcell,drop=F]
      percellanno <- cellanno[sampcell,,drop=F]
      perpsn <- perpsn[sampcell]
      colnames(perexpr) <- percellanno[,1] <- names(perpsn) <- paste0('cell_',1:length(perpsn))
      ### runalgo 
      perfit <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=perpsn, design=design, knotnum=NULL, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, parallel=F)
      perfit <- sapply(perfit$parameter,function(i) unname(i$ll))
      perfit <- perfit[row.names(expr)]
    })
  } else {
    perll <- mclapply(1:permuiter,function(did) {
      ### permute time within each sample , 20200630
      perpsn <- sapply(rownames(design), function(s){
        tmpid <- cellanno[cellanno[,2] == s, 1]  # subset cells
        tmppsn <- pseudotime[names(pseudotime) %in% tmpid] # subset time
        names(tmppsn) <- sample(names(tmppsn)) # permute time
        tmppsn
      })
      names(perpsn) <- NULL
      perpsn <- unlist(perpsn)
      perpsn <- perpsn[colnames(expr)]
      ### bootstrap within each sample , 20200630
      sampcell <- as.vector(unlist(lapply(unique(cellanno[,2]), function(p){
        sample(which(cellanno[,2] == p), replace = T)
      })))
      perexpr <- expr[,sampcell,drop=F]
      percellanno <- cellanno[sampcell,,drop=F]
      perpsn <- perpsn[sampcell]
      colnames(perexpr) <- percellanno[,1] <- names(perpsn) <- paste0('cell_',1:length(perpsn))
      ### run algo
      perfit <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=prepsn, design=design, knotnum=NULL, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, parallel = F)
      perfit <- sapply(perfit$parameter,function(i) unname(i$ll))
      perfit <- perfit[row.names(expr)]
    }, mc.cores = ncores)
  }
  perll <- do.call(cbind,perll)
  pval <- sapply(1:nrow(perll), function(i) pnorm(orill[i],mean(perll[i,]),sd(perll[i,]),lower.tail = F))
  fdr <- p.adjust(pval,method='fdr')
  names(fdr) <- row.names(perll)
  foldchange <- orill - rowMeans(perll)
  return(list(fdr = fdr, perll = perll, knotnum = knotnum, foldchange = foldchange, parameter=orifit$parameter, pseudotime = pseudotime))
}

