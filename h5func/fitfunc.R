fitfunc <- function(iter, diffType = 'overall', gene = NULL, test.type = 'Time', EMmaxiter=100, EMitercutoff=0.1, verbose=F, ncores=1, expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design) {
  ## this function serves the function testpt().
  ## return(list(fitres.full = fitres.full, fitres.null = fitres.null))
  ## ncores = 1 or otherwise meaningless since the upper function is running in parallel
  ## expr: hdf5 path and file name
  print(paste0('iter ', iter, '\n'))
  if (toupper(test.type)=='TIME') {
    if (iter == 1){
      fitres.full <- fitpt(expr=expr, pseudotime=pseudotime, design=design[,1,drop=FALSE],targetgene=gene, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores, model=1)
      fitres.null <- fitpt.m0(expr=expr, pseudotime=pseudotime, design=design[,1,drop=FALSE],targetgene=gene, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose)
      return(list(fitres.full = fitres.full, fitres.null = fitres.null))
    } else {
      perpsn <- lapply(rownames(design), function(s){
        tmpid <- cellanno[cellanno[,2] == s, 1]  # subset cells
        tmppsn <- pseudotime[names(pseudotime) %in% tmpid] # subset time
        names(tmppsn) <- sample(names(tmppsn)) # permute time
        tmppsn
      })
      names(perpsn) <- NULL
      perpsn <- unlist(perpsn)
      sampcell <- as.vector(unlist(lapply(unique(cellanno[,2]), function(p){
        sample(which(cellanno[,2] == p), replace = T)
      }))) ### bootstrap
      percellanno <- cellanno[sampcell,,drop=F]
      perpsn <- perpsn[names(pseudotime)]
      perpsn <- perpsn[sampcell]
      boot <- data.frame(percellanno[,1],paste0('cell_',1:length(perpsn)),stringsAsFactors = F) #### rename cells
      percellanno[,1] <- names(perpsn) <- paste0('cell_',1:length(perpsn)) ## save the original cell name and permuted cell names relation, not used here actually, because hdf5 file already save the cells for each sample seperately
      ####
      tryCatch(fitres.full <- fitpt(expr=expr, pseudotime=perpsn, design=design, boot=boot,targetgene=gene,EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model = 1), warning = function(w){}, error = function(e) {})
      tryCatch(fitres.null <- fitpt.m0(expr=expr, pseudotime=perpsn, design=design, boot=boot,targetgene=gene,EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose), warning = function(w){}, error = function(e) {})
      if (exists('fitres.full') & exists('fitres.null')) {
        print(paste0('iter ', iter, ' success!'))
        return(list(fitres.full = fitres.full, fitres.null = fitres.null))
      } else {
        print(paste0('iter ', iter, ' try again!'))
        return(NULL)
      }
    }
  } else if (toupper(test.type)=='VARIABLE'){
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
      fitres.full <- fitpt(expr,  pseudotime, design,targetgene=gene, maxknotallowed=10, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model = mod.full)
      fitres.null <- fitpt(expr,  pseudotime, design,targetgene=gene, maxknotallowed=10, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=1, model = mod.null, knotnum = fitres.full[[2]])
      if (exists('fitres.full') & exists('fitres.null')) {
        print(paste0('iter ', iter, ' success!'))
        return(list(fitres.full = fitres.full, fitres.null = fitres.null))
      } else if (!exists('fitres.full')){
        print(paste0('iter ', iter, 'full model failed. Skip ...'))
        return(NULL)
      } else {
        print(paste0('iter ', iter, 'null model failed. Skip ...'))
        return(NULL)
      }
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
    sampcell <- sample(1:length(pseudotime),replace=T) ## boostrap cells
    percellanno <- cellanno[sampcell,,drop=F]
    psn <- pseudotime
    psn <- psn[sampcell]
    boot <- data.frame(percellanno[,1],paste0('cell_',1:length(psn)),stringsAsFactors = F)
    percellanno[,1] <- names(psn) <- paste0('cell_',1:length(psn))
    
    
    fitres.full <- fitpt(expr,  psn, perdesign, boot=boot,targetgene=gene,maxknotallowed=10, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores, model = mod.full)
    fitres.null <- fitpt(expr,  psn, perdesign, boot=boot,targetgene=gene,maxknotallowed=10, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, verbose=verbose, ncores=ncores, model = mod.null, knotnum = fitres.full[[2]])
    if (exists('fitres.full') & exists('fitres.null')) {
      print(paste0('iter ', iter, ' success!'))
      return(list(fitres.full = fitres.full, fitres.null = fitres.null))
    } else if (!exists('fitres.full')){
      print(paste0('iter ', iter, 'full model failed. Skip ...'))
      return(NULL)
    } else {
      print(paste0('iter ', iter, 'null model failed. Skip ...'))
      return(NULL)
    }
  }
}
}



