testpt <- function(expr, cellanno, pseudotime, design=NULL, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type='Time', test.pattern = 'overall', test.position = 'all', fit.resolution = 1000, return.all.data = TRUE, demean = TRUE) {
  set.seed(12345)
  library(splines)
  cellanno = data.frame(Cell = as.character(cellanno[,1]), Sample = as.character(cellanno[,2]), stringsAsFactors = FALSE)
  expr.ori <- expr[, names(pseudotime)]
  cellanno <- cellanno[match(names(pseudotime), cellanno[,1]), ]
  if (demean){
    ## demean
    expr.demean <- lapply(unique(cellanno[,2]), function(s){
      tmp <- expr.ori[, cellanno[cellanno[,2] == s, 1]]
      tmp2 <- tmp- rowMeans(tmp)
    })
    expr.demean <- do.call(cbind, expr.demean)
    expr <- expr.demean[, colnames(expr.ori)]
  } else {
    expr <- expr.ori
  }
  if (type=='Time') {
    unis <- unique(cellanno[,2])
    design = matrix(1,nrow=length(unis),ncol=1,dimnames = list(unis,'intercept'))
  } else {
    design = as.matrix(design)  
  }
  fitfunc <- function(iter) {
    if (iter==1) {
      fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, test.pattern = test.pattern, test.position = test.position)  
    } else {
      if (type=='Time') {
        print('testing Time ...')
        print(paste0('iteration ', iter))
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
        # << --- 20200926 
        # fitpt(expr=perexpr, cellanno=percellanno, pseudotime=perpsn, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, test.pattern = test.pattern, test.position = test.position)
        tryCatch(res <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=perpsn, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, test.pattern = test.pattern, test.position = test.position), warning = function(w){}, error = function(e) {})
        if (exists('res')) {
          return(res)
        } else {
          return(NULL)
        }
        # --- 20200926 --->>
      } else if (type=='Variable'){
        print('testing Variable ...')
        print(paste0('iteration ', iter))
        dn <- paste0(as.vector(design),collapse = '_')
        perdn <- dn
        while(perdn==dn) {
          perid <- sample(1:nrow(design))
          perdesign <- design[perid,,drop=F]
          perdn <- paste0(as.vector(perdesign),collapse = '_')  
        }
        row.names(perdesign) <- row.names(design)
        sampcell <- sample(1:ncol(expr),replace=T) ## boostrap cells
        perexpr <- expr[,sampcell,drop=F]
        percellanno <- cellanno[sampcell,,drop=F]
        psn <- pseudotime[colnames(expr)]
        psn <- psn[sampcell]
        colnames(perexpr) <- percellanno[,1] <- names(psn) <- paste0('cell_',1:length(psn))
        fitpt(expr=perexpr, cellanno=percellanno, pseudotime=psn, design=perdesign, ori.design = design, test.pattern = test.pattern, test.position = test.position, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1)  
      }
    }
  }
  print('fitting model ...')
  if (ncores == 1){
    fit <- lapply(1:(permuiter+1),fitfunc)
  } else {
    fit <- mclapply(1:(permuiter+1),function(i){set.seed(i); fitfunc(i)}, mc.cores = ncores)
  }
  # << --- 20200926 
  permuiter <- sum(!sapply(fit,is.null)) -1
  fit <- fit[!sapply(fit,is.null)]
  # --- 20200926 --- >>
  orifit <- fit[[1]]
  knotnum <- orifit$knotnum
  orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)[row.names(expr)]
  print('performing permutation test ...')
  perll <- sapply(2:(permuiter+1),function(i) sapply(fit[[i]]$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
  print('calculationg p-values ...')
  pval <- sapply(1:nrow(perll), function(i) {
    z <- perll[i,]
    den <- density(z,bw='SJ')$bw
    mean(pnorm(orill[i], z, sd=den,lower.tail = F))
  })
  fdr <- p.adjust(pval,method='fdr')
  names(pval) <- names(fdr) <- row.names(perll)
  foldchange <- orill - rowMeans(perll)
  #---------------------------------
  # use beta2 to get mean difference
  # --------------------------------
  max.abs.beta2 <- sapply(names(orifit$parameter), function(g){
      beta <- orifit$parameter[[g]]$beta
      tmp <- abs(unlist(sapply(1:length(beta), function(i){
        if (i%%2 == 0) beta[i]
      })))
      a <- ceiling(length(tmp)/3)
      b <- ceiling(length(tmp)*2/3)
      if (test.pattern == 'overall'){
        if (test.position == 'all'){
          return(max(tmp, na.rm = TRUE))
        } else if (test.position == 'start'){
          return(max(tmp[c(1, seq(2, a))], na.rm = TRUE))
        } else if (test.position == 'middle'){
          return(max(tmp[c(1, seq(a+1, b))], na.rm = TRUE))
        } else if (test.position == 'end'){
          return(max(tmp[c(1, seq(b+1, length(tmp)))], na.rm = TRUE))
        }
      } else if (test.pattern == 'slope'){
        if (test.position == 'all'){
          return(max(tmp[seq(2, length(tmp))], na.rm = TRUE))
        } else if (test.position == 'start'){
          return(max(tmp[seq(2, a)], na.rm = TRUE))
        } else if (test.position == 'middle'){
          return(max(tmp[seq(a+1, b)], na.rm = TRUE))
        } else if (test.position == 'end'){
          return(max(tmp[seq(b+1, length(tmp))], na.rm = TRUE))
        }
      } else if (test.pattern == 'intercept'){
        return(tmp[1])
      }
  })
  # -------------------------------
  pred <- predict_fitting(expr = expr,knotnum = knotnum, design = design, cellanno = cellanno, pseudotime = pseudotime[colnames(expr)])
  if (return.all.data){
    if (demean){
      return(list(fdr = fdr, foldchange = foldchange, pvalue = pval, max.abs.beta2 = max.abs.beta2, parameter=orifit$parameter, orill=orill, perll = perll, knotnum = knotnum,  pseudotime = pseudotime[colnames(expr)], predict.values = pred[,colnames(expr)], design = design, cellanno = cellanno, expr.demean = expr, expr.ori = expr.ori))
    } else {
      return(list(fdr = fdr, foldchange = foldchange, pvalue = pval, max.abs.beta2 = max.abs.beta2, parameter=orifit$parameter, orill=orill, perll = perll, knotnum = knotnum,  pseudotime = pseudotime[colnames(expr)], predict.values = pred[,colnames(expr)], design = design, cellanno = cellanno, expr.ori = expr))
    } 
  } else {
    return(list(fdr = fdr, foldchange = foldchange, pvalue = pval, max.abs.beta2 = max.abs.beta2, parameter=orifit$parameter, orill=orill, perll = perll, knotnum = knotnum, predict.values = pred[,colnames(expr)]))
  } 
}






