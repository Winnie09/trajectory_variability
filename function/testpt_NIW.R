testpt <- function(expr, cellanno, pseudotime, design=NULL, permuiter=100, EMmaxiter=1000, EMitercutoff=0.01, verbose=F, ncores=detectCores(), type='Time', test.pattern = 'overall', test.position = 'all', fit.resolution = 1000, return.all.data = TRUE, demean = FALSE) {
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
  fitfunc <- function(iter, diffType = 'both') {
    print(paste0('iter ', iter, '\n'))
    if (iter==1) {
      if (diffType == 'both' | diffType == 'trendDiff') {
        model = 3
      } else if (diffType == 'meanDiff') {
        model = 2
      } 
      fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, test.pattern = 'overall', test.position = test.position, model = model) 
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
        if (diffType == 'both'){
          test.pattern = 'overall'
          model = 3
        } else if (diffType == 'meanDiff'){
          test.pattern = 'intercept'
          model = 2
        } else if (diffType == 'trendDiff'){
          test.pattern = 'slope'
          model = 3
        }
        fitpt(expr=perexpr, cellanno=percellanno, pseudotime=psn, design=perdesign, ori.design = design, test.pattern = test.pattern, test.position = test.position, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model = model)  
      }
    }
  }
  print('fitting model ...')
  if (type == 'Variable'){
    ## meanDiff pvalues: Model 2 vs. model 1
    if (ncores == 1){
      fit <- lapply(1:(permuiter+1),function(i) fitfunc(iter = i, diffType = 'meanDiff'))
    } else {
      fit <- mclapply(1:(permuiter+1),function(i){set.seed(i); fitfunc(iter = i, diffType = 'meanDiff')}, mc.cores = ncores)
    }
    permuiter <- sum(!sapply(fit,is.null)) -1
    fit <- fit[!sapply(fit,is.null)]
    orifit <- fit[[1]]
    knotnum <- orifit$knotnum
    orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)[row.names(expr)]
    print('meandiff: performing permutation test ...')
    perll <- sapply(2:(permuiter+1),function(i) sapply(fit[[i]]$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
    print('meanDiff: calculationg p-values ...')
    pval <- sapply(1:nrow(perll), function(i) {
      z <- perll[i,]
      den <- density(z,bw='SJ')$bw
      mean(pnorm(orill[i], z, sd=den,lower.tail = F))
    })
    fdr <- p.adjust(pval,method='fdr')
    names(pval) <- names(fdr) <- row.names(perll)
    foldchange <- orill - rowMeans(perll)
    res <- data.frame(meanDiff.fdr = fdr, meanDiff.fc = foldchange, meanDiff.pvalue = pval, stringsAsFactors = FALSE)
    
    ## trendDiff pvalues: Model 3 vs. model 2
    if (ncores == 1){
      fit <- lapply(1:(permuiter+1), function(i) fitfunc(iter = i, diffType = 'trendDiff'))
    } else {
      fit <- mclapply(1:(permuiter+1), function(i){set.seed(i); fitfunc(iter = i, diffType = 'trendDiff')}, mc.cores = ncores)
    }
    permuiter <- sum(!sapply(fit,is.null)) -1
    fit <- fit[!sapply(fit,is.null)]
    orifit <- fit[[1]]
    knotnum <- orifit$knotnum
    orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)[row.names(expr)]
    print('trendDiff: performing permutation test ...')
    perll <- sapply(2:(permuiter+1),function(i) sapply(fit[[i]]$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
    print('trendDiff: calculationg p-values ...')
    pval <- sapply(1:nrow(perll), function(i) {
      z <- perll[i,]
      den <- density(z,bw='SJ')$bw
      mean(pnorm(orill[i], z, sd=den,lower.tail = F))
    })
    fdr <- p.adjust(pval,method='fdr')
    names(pval) <- names(fdr) <- row.names(perll)
    foldchange <- orill - rowMeans(perll)
    res <- cbind(res, trendDiff.fdr = fdr, trendDiff.fc = foldchange, trendDiff.pvalue = pval)
  } 
    

  ## overall (both) pvalues: Model 3 vs. model 2
  if (ncores == 1){
    fit <- lapply(1:(permuiter+1),function(i) fitfunc(iter = i, diffType = 'both'))
  } else {
    fit <- mclapply(1:(permuiter+1),function(i){set.seed(i); fitfunc(iter = i, diffType = 'both')}, mc.cores = ncores)
  }
  permuiter <- sum(!sapply(fit,is.null)) -1
  fit <- fit[!sapply(fit,is.null)]
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
  if (type == 'Variable'){
    res <- cbind(res, both.fdr = fdr, both.fc = foldchange, both.pvalue = pval)
  } else if (type == 'Time'){
    res <- data.frame(fdr = fdr, fc = foldchange, pvalue = pval, stringsAsFactors = FALSE)
  }
    
  pred <- predict_fitting(expr = expr,knotnum = knotnum, design = design, cellanno = cellanno, pseudotime = pseudotime[colnames(expr)])
  if (return.all.data){
    if (demean){
      return(list(statistics = res, parameter=orifit$parameter, orill=orill, perll = perll, knotnum = knotnum,  pseudotime = pseudotime[colnames(expr)], predict.values = pred[,colnames(expr)], design = design, cellanno = cellanno, expr.demean = expr, expr.ori = expr.ori))
    } else {
      return(list(statistics = res, parameter=orifit$parameter, orill=orill, perll = perll, knotnum = knotnum,  pseudotime = pseudotime[colnames(expr)], predict.values = pred[,colnames(expr)], design = design, cellanno = cellanno, expr.ori = expr))
    } 
  } else {
    return(list(statistics = res, parameter=orifit$parameter, orill=orill, perll = perll, knotnum = knotnum, predict.values = pred[,colnames(expr)]))
  } 
}





