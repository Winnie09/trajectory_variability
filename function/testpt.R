testpt <- function(expr, cellanno, pseudotime, design=NULL, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Time') {
  if (is.data.frame(pseudotime)) {
    pseudotime <- data.frame(Cell = pseudotime[,1], Pseudotime = as.numeric(pseudotime[,2]), stringsAsFactors = FALSE)
    pseudotime <- pseudotime[order(pseudotime[,2]), ]
    psn <- pseudotime[,2]
    names(psn) <- pseudotime[,1]
    pseudotime <- psn
  }
  set.seed(12345)
  expr <- expr[, names(pseudotime)]
  cellanno <- cellanno[match(names(pseudotime), cellanno[,1]), ]
  if (type=='Time') {
    unis <- unique(cellanno[,2])
    design = matrix(1,nrow=length(unis),ncol=1,dimnames = list(unis,NULL))
  }
  
  fitfunc <- function(iter) {
    if (iter==1) {
      fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)  
    } else {
      if (type=='Time') {
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
        fitpt(expr=perexpr, cellanno=percellanno, pseudotime=perpsn, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)
      } else if (type=='Variable'){
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
    den <- density(z,bw='SJ')$bw
    mean(pnorm(orill[i], z, sd=den,lower.tail = F))
  })
  fdr <- p.adjust(pval,method='fdr')
  names(fdr) <- row.names(perll)
  foldchange <- orill - rowMeans(perll)
  return(list(fdr = fdr, orill=orill, perll = perll, knotnum = knotnum, foldchange = foldchange, parameter=orifit$parameter, pseudotime = pseudotime))
}


