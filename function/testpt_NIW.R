# permuiter=100; EMmaxiter=1000; EMitercutoff=0.01; verbose=F; ncores=detectCores(); test.type='Variable'; fit.resolution = 1000; return.all.data = TRUE; demean = FALSE; overall.only = T; test.method = 'EM'

testpt <- function(expr, cellanno, pseudotime, design=NULL, permuiter=100, EMmaxiter=1000, EMitercutoff=0.01, verbose=F, ncores=detectCores(), test.type='Time', fit.resolution = 1000, return.all.data = TRUE, demean = FALSE, overall.only = F, test.method = 'EM') {
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
  if (test.type=='Time') {
    unis <- unique(cellanno[,2])
    design = matrix(1,nrow=length(unis),ncol=1,dimnames = list(unis,'intercept'))
  } else {
    design = as.matrix(design)  
  }
  
  if (test.method == 'chisq' & test.type == 'Variable'){
    res3 <- fitpt(expr, cellanno, pseudotime, design, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.1, verbose=F, ncores=1, model = 3)
    ll3 <- sapply(res3$parameter,function(i) i$ll)
    res2 <- fitpt(expr, cellanno, pseudotime, design, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.1, verbose=F, ncores=1, model = 3)
    ll2 <- sapply(res2$parameter,function(i) i$ll)
    res1 <- fitpt(expr, cellanno, pseudotime, design=design, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.1, verbose=F, ncores=1, model = 1)
    ll1 <- sapply(res1$parameter,function(i) i$ll)
    
    paradiff31 <- sapply(res3$parameter,function(i) length(i$beta))-sapply(res1$parameter,function(i) length(i$beta))
    paradiff32 <- sapply(res3$parameter,function(i) length(i$beta))-sapply(res2$parameter,function(i) length(i$beta))
    paradiff21 <- sapply(res2$parameter,function(i) length(i$beta))-sapply(res1$parameter,function(i) length(i$beta))
    lldiff <- ll3-ll1
    
    pval.chisq.overall <- pchisq(2*(ll3-ll1),df=paradiff31,lower.tail = F)
    fdr.chisq.overall <- p.adjust(pval.chisq.overall, method='fdr')
    pval.chisq.trendDiff <- pchisq(2*(ll3-ll2),df=paradiff32,lower.tail = F)
    fdr.chisq.trendDiff <- p.adjust(pval.chisq.trendDiff, method='fdr')
    pval.chisq.meanDiff <- pchisq(2*(ll2-ll1),df=paradiff21,lower.tail = F)
    fdr.chisq.meanDiff <- p.adjust(pval.chisq.meanDiff, method='fdr')
    
    res <- data.frame(fdr.chisq.overall = fdr.chisq.overall, pval.chisq.overall = pval.chisq.overall,
                      fdr.chisq.trendDiff = fdr.chisq.trendDiff, pval.chisq.trendDiff = pval.chisq.trendDiff, 
                      fdr.chisq.meanDiff = fdr.chisq.meanDiff, pval.chisq.meanDiff = pval.chisq.meanDiff,
                      stringsAsFactors = FALSE)
    return(list(statistics = res, ll1 = ll1, ll2 = ll2, ll3 = ll3))  ## function return
  } else if (test.method == 'EM'){
    fitfunc <- function(iter, diffType = 'overall', gene = rownames(expr)) {
      expr <- expr[gene, ,drop=FALSE]
      print(paste0('iter ', iter, '\n'))
        if (test.type=='Time') {
          if (iter == 1){
            fitres.full <- fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model=1)
            fitres.null <- fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model=0)
            return(list(fitres.full = fitres.full, fitres.null = fitres.null))
          } else {
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
            tryCatch(fitres.full <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=psn, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model = 1), warning = function(w){}, error = function(e) {})
            tryCatch(fitres.null <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=psn, design=design, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model = 0), warning = function(w){}, error = function(e) {})
            if (exists('fitres.full') & exists('fitres.null')) {
              return(list(fitres.full = fitres.full, fitres.null = fitres.null))
            } else {
              return(NULL)
            }
          }
        } else if (test.type=='Variable'){
          print('testing Variable ...')
          if (diffType == 'overall'){
            mod.full = 3
            mod.null = 1
          } else if (diffType == 'meanDiff'){
            mod.full = 2
            mod.null = 1
          } else if (diffType == 'trendDiff'){
            mod.full = 3
            mod.null = 2
          }
          if (iter == 1){
             fitres.full <- fitpt(expr, cellanno, pseudotime, design, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.1, verbose=F, ncores=1, model = mod.full)
             fitres.null <- fitpt(expr, cellanno, pseudotime, design, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.1, verbose=F, ncores=1, model = mod.null)
             if (exists('fitres.full') & exists('fitres.null')) {
                return(list(fitres.full = fitres.full, fitres.null = fitres.null))
              } else {
                return(NULL)
              }
          } else {
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
            fitres.full <- fitpt(perexpr, percellanno, psn, perdesign, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.1, verbose=F, ncores=1, model = mod.full)
            fitres.null <- fitpt(perexpr, percellanno, psn, perdesign, maxknotallowed=10, EMmaxiter=1000, EMitercutoff=0.1, verbose=F, ncores=1, model = mod.null)
           if (exists('fitres.full') & exists('fitres.null')) {
              return(list(fitres.full = fitres.full, fitres.null = fitres.null))
            } else {
              return(NULL)
            }
          }
        }
    }
    
    print('fitting model ...')
    ## overall (overall) pvalues: Model 3 vs. model 1
    if (ncores == 1){
      fit <- lapply(1:(permuiter+1),function(i) fitfunc(iter = i, diffType = 'overall'))
    } else {
      fit <- mclapply(1:(permuiter+1),function(i){set.seed(i); fitfunc(iter = i, diffType = 'overall')}, mc.cores = ncores)
    }
    fit <- fit[!sapply(fit,is.null)]
    knotnum <- fit[[1]]$fitres.full$knotnum
    parameter <- fit[[1]]$fitres.full$parameter
    ll.full <- sapply(1:(length(fit)),function(i) sapply(fit[[i]]$fitres.full$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
    ll.null <- sapply(1:(length(fit)),function(i) sapply(fit[[i]]$fitres.null$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
    llr.overall <- ll.full - ll.null
    pval.overall <- sapply(1:nrow(llr.overall), function(i) {
      z <- llr.overall[i,2:ncol(llr.overall)]
      den <- density(z)$bw
      mean(pnorm(llr.overall[i,1], z, sd=den, lower.tail = F))
    })
    fdr.overall <- p.adjust(pval.overall,method='fdr')
    names(pval.overall) <- names(fdr.overall) <- row.names(llr.overall)
    z.score <- (llr.overall[,1] - rowMeans(llr.overall[,2:(ncol(llr.overall))]))/apply(llr.overall[,2:(ncol(llr.overall))],1,sd)
    res.overall <- data.frame(fdr.overall = fdr.overall, z.overall = z.score, pvalue.overall = pval.overall, stringsAsFactors = FALSE)
    
    if (test.type == 'Variable' & !overall.only){
      ## meanDiff pvalues: Model 2 vs. model 1
      if (ncores == 1){
        fit <- lapply(1:(permuiter+1),function(i) fitfunc(iter = i, diffType = 'meanDiff', gene = names(fdr.overall)[fdr.overall<0.05]))
      } else {
        fit <- mclapply(1:(permuiter+1),function(i){set.seed(i); fitfunc(iter = i, diffType = 'meanDiff', gene = names(fdr.overall)[fdr.overall<0.05])}, mc.cores = ncores)
      }
      fit <- fit[!sapply(fit,is.null)]
      ll.full <- sapply(1:(length(fit)),function(i) sapply(fit[[i]]$fitres.full$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      ll.null <- sapply(1:(length(fit)),function(i) sapply(fit[[i]]$fitres.null$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      llr <- ll.full - ll.null
      if (sum(fdr.overall<0.05) == 1){
        z <- llr[2:length(llr)]
        den <- density(z)$bw
        fdr <- pval <- mean(pnorm(llr[1], z, sd=den,lower.tail = F))
        names(pval) <- names(fdr) <- names(fdr.overall)[fdr.overall<0.05]
        z.score <- (llr[1] - mean(z))/sd(z)
      } else {
        pval <- sapply(1:nrow(llr), function(i) {
          z <- llr[i,2:ncol(llr)]
          den <- density(z)$bw
          mean(pnorm(llr[i,1], z, sd=den,lower.tail = F))
        })
        fdr <- p.adjust(pval,method='fdr')
        names(pval) <- names(fdr) <- row.names(perll)
        z.score <- (llr[,1] - rowMeans(llr[,2:(ncol(llr))]))/apply(llr[,2:(ncol(llr))],1,sd)
      }
      res.meanDiff <- data.frame(fdr.meanDiff = fdr, z.meanDiff = z.score, pvalue.meanDiff = pval, stringsAsFactors = FALSE)
      
      ## trendDiff pvalues: Model 3 vs. model 2
      if (ncores == 1){
        fit <- lapply(1:(permuiter+1), function(i) fitfunc(iter = i, diffType = 'trendDiff', gene = names(fdr.overall)[fdr.overall<0.05]))
      } else {
        fit <- mclapply(1:(permuiter+1), function(i){set.seed(i); fitfunc(iter = i, diffType = 'trendDiff', gene = names(fdr.overall)[fdr.overall<0.05])}, mc.cores = ncores)
      }
      fit <- fit[!sapply(fit,is.null)]
      ll.full <- sapply(1:(length(fit)),function(i) sapply(fit[[i]]$fitres.full$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      ll.null <- sapply(1:(length(fit)),function(i) sapply(fit[[i]]$fitres.null$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      llr <- ll.full - ll.null
      if (sum(fdr.overall<0.05) == 1){
        z <- llr[2:length(llr)]
        den <- density(z)$bw
        fdr <- pval <- mean(pnorm(llr[1], z, sd=den,lower.tail = F))
        names(pval) <- names(fdr) <- names(fdr.overall)[fdr.overall<0.05]
        z.score <- (llr[1] - mean(z))/sd(z)
      } else {
        pval <- sapply(1:nrow(llr), function(i) {
          z <- llr[i,2:ncol(llr)]
          den <- density(z)$bw
          mean(pnorm(llr[i,1], z, sd=den,lower.tail = F))
        })
        fdr <- p.adjust(pval,method='fdr')
        names(pval) <- names(fdr) <- row.names(perll)
        z.score <- (llr[,1] - rowMeans(llr[,2:(ncol(llr))]))/apply(llr[,2:(ncol(llr))],1,sd)
      }
      res.trendDiff <- data.frame(fdr.trendDiff = fdr, z.trendDiff = z.score, pvalue.trendDiff = pval, stringsAsFactors = FALSE)
      res <- matrix(NA, nrow = nrow(res.overall), ncol = 9, 
                      dimnames = list(rownames(res.overall),c(colnames(res.overall), colnames(res.meanDiff), colnames(res.trendDiff))))
      res[rownames(res.overall), colnames(res.overall)] <- as.matrix(res.overall)
      res[rownames(res.trendDiff), colnames(res.trendDiff)] <- as.matrix(res.trendDiff)
      res[rownames(res.meanDiff), colnames(res.meanDiff)] <- as.matrix(res.meanDiff)
    } else if (test.type == 'Variable' & overall.only){
      res <- res.overall
    }
      
    return(list(statistics = res, 
                parameter=parameter, 
                ll.full = ll.full,
                ll.null = ll.null,
                knotnum = knotnum))
    
  }
}



