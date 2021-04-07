cellPropTest <- function(cellanno, pseudotime, design=NULL, ncores=detectCores()) {
  ptw <- cut(pseudotime,seq(min(pseudotime),max(pseudotime),length.out = 100),include.lowest = T)
  ptdat <- table(ptw,cellanno[match(names(pseudotime),cellanno[,1]),2])
  ptdat <- t(t(ptdat)/colSums(ptdat)) ## divided by rowsum (rowsum = 1). interval * samples. 
  ptdat <- as.data.frame(ptdat)
  colnames(ptdat) <- c('pt','s','prop')
  ptdat[,1] <- match(ptdat[,1],levels(ptw))
  
  ptdat$cell <- paste0('cell',1:nrow(ptdat))
  ptexpr <- t(ptdat[,c('prop','prop'),drop=F])
  colnames(ptexpr) <- ptdat$cell
  
  ptpt <- ptdat$pt
  names(ptpt) <- ptdat$cell
  res <- testpt(expr=ptexpr, cellanno=data.frame(cell=ptdat$cell,sample=ptdat$s), pseudotime=ptpt, design=design, ncores=ncores, test.type = 'Variable', demean = FALSE, test.method = 'permutation', ncores.fit = 1,fix.all.zero=F)
  res$statistics <- res$statistics[1, 2:3]
  return(res)
}

