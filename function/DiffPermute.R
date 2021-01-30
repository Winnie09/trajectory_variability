DiffPermute <- function(GeneByCellExpr, cellanno, pseudotime, design, parallel=FALSE, maxCore = 4, maxCell = NULL){
  if (!is.null(maxCell)){
    selectcell <- lapply(unique(cellanno[,2]), function(s){
      tmp <- cellanno[cellanno[,2]==s, , drop=F]
      if (nrow(tmp) > maxCell){
        set.seed(12345)
        tmp[sort(sample(seq(1,nrow(tmp)), maxCell)), ,drop=F]
      } else {
        tmp
      }
        
    })
    cellanno <- do.call(rbind, selectcell)
  }
  GeneByCellExpr <- GeneByCellExpr[, cellanno[,1]]
  design = cbind(1,design)
  library(gtools)
  # per <- permutations(nrow(design),nrow(design)) ## use when small #samples 
  per <- t(sapply(1:100, function(i){
    set.seed(i)
    sample(rownames(design))
  }))
  per <- rbind(per,rownames(design)) ## in case it is a continuours and hard to get repetitive permuted result
  per <- unique(per)  ## in case one permutation is the same as the design
  id <- which(!duplicated(apply(per,1,function(i) paste0(as.vector(design[i,]),collapse = '_'))))
  per <- per[id,]
  oriid <- which(apply(per,1,function(i) paste0(as.vector(design[i,]),collapse = '_'))==paste0(as.vector(design),collapse = '_'))
  library(parallel)
  if (parallel){
    perll <- mclapply(1:nrow(per),function(i) {
      print(i)
      perdesign <- design[per[i,],,drop=F]
      row.names(perdesign) <- row.names(design)
      diffpt(expr=GeneByCellExpr,design=perdesign,pseudotime=pseudotime,num.knot = 3,cellanno = cellanno, verbose = T)$logL
    }, mc.cores = maxCore)
  } else {
    perll <- lapply(1:nrow(per),function(i) {
      print(i)
      perdesign <- design[per[i,],,drop=F]
      row.names(perdesign) <- row.names(design)
      diffpt(expr=GeneByCellExpr,design=perdesign,pseudotime=pseudotime,num.knot = 3,cellanno = cellanno, verbose = T)$logL
    })
  } 
  perll <- do.call(cbind,perll)
  pval <- sapply(1:nrow(perll), function(i) pnorm(perll[i,oriid],mean(perll[i,-oriid]),sd(perll[i,-oriid]),lower.tail = F))
  fdr <- p.adjust(pval,method='fdr')
  names(fdr) <- row.names(perll)
  res <- data.frame(P.Value = pval, adj.P.Val = fdr, stringsAsFactors = F)
  rownames(res) <- names(fdr)
  res <- res[order(res[,2]),]
  final <- list()
  final[['res']] <- res
  final[['perll']] <- perll
  return(final)
}
