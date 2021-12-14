saveh5 <- function(expr,pseudotime,cellanno,path,fix.all.zero = TRUE, cutoff = 1e-5, sd.adjust = 1) {
## path: file location of the  hdf5 file
## 
  expr <- expr[,names(pseudotime)]
  if (fix.all.zero){
    sdm <- sapply(unique(cellanno[,2]),function(us) {
      tmp <- expr[,cellanno[,2]==us, drop=FALSE]
      m <- rowMeans(tmp)
      rowMeans(tmp*tmp)-m*m
    })<cutoff ## first version <20210407 is 0.
    gid <- which(rowSums(sdm) > 0)  ## identify if any genes have sd=0 expression in any one of the samples
    if (length(gid) > 0) { ## if yes, for those genes, add a white-noise with sd=1e-5 on the sample with sd=0.
      mask <- sdm[gid,rep(1:ncol(sdm),as.vector(table(cellanno[,2])[colnames(sdm)])),drop=F]
      colnames(mask) <- unlist(sapply(colnames(sdm),function(i) cellanno[cellanno[,2]==i,1]))
      expr[gid,] <- expr[gid,] + mask[,colnames(expr),drop=F] * matrix(rnorm(length(mask),sd=sd.adjust), nrow=length(gid)) ## before 20211126, sd.adjust = 1e-5 
      rm('mask')
    }
  }
  
  h5createFile(path) ## create a hdf5 file
  for (s in unique(cellanno[,2])) { ## operation for each sample
    h5createGroup(path,s) ## create a group for each sample
    d <- expr[,cellanno[cellanno[,2]==s,1]] 
    h5createDataset(path, paste0(s,"/expr"), c(nrow(d),ncol(d)), chunk=c(1,ncol(d)), level=1) ## demonstrate will create a dataset of size nrow*ncol, chunk is the smallest saving unit
    h5write(d, file=path, paste0(s,"/expr"), index=list(1:nrow(d),NULL))  ## write the xpr. index: allow to get the data for each gene across all cell (1:nrow(d)), not allow to get a cell across all genes (second arg of index is NULL). The matrix saved here has no row or col names.
    
    h5createDataset(path, paste0(s,"/feature"), c(nrow(d)), chunk=c(nrow(d)), level=1,storage.mode='character',size=max(nchar(rownames(d))))
    h5write(rownames(d), file=path, name=paste0(s,"/feature"), index=list(1:nrow(d))) ## save feature names (row names, gene names)
    
    h5createDataset(path, paste0(s,"/barcode"), c(ncol(d)), chunk=c(ncol(d)), level=1,storage.mode='character',size=max(nchar(colnames(d))))
    h5write(colnames(d), file=path, name=paste0(s,"/barcode"), index=list(1:ncol(d))) ## save barcode names (col names, cell names).
  }
  H5close()
}


