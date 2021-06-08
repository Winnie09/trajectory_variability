cluster_gene <- function(testobj, 
                         gene,
                         k,
                         k.auto = FALSE,
                         type = 'time',
                         method = 'kmeans', 
                         scale.difference = F){
  ## testobj: output object from testpt().
  ## gene: a character vector of genes.
  ## k: number of clusters. 2 are for meanSig, and the remaining k-2 clusters are for DDG of NOT meanSig.
  ## type: "time" or "variable". A character denoting the population fit is on "time" or "variable".  Default is "time". ## method: 'kmeans'(default) or 'hierarchical'.
  ## scale.difference: if FALSE (default), do not standardize the group difference. If TRUE, scale the group difference by the maximum of absolute values. 
  if (type == 'time'){
    if ('populationFit' %in% names(testobj)) {
      fit <- testobj$populationFit
    } else {
      fit <- getPopulationFit(testobj, gene = gene, type = type)
    }
  } else if (type == 'variable'){
    if ('covariateGroupDiff' %in% names(testobj)){
      fit <- testobj$covariateGroupDiff
    } else{
      fit <- getCovariateGroupDiff(testobj = testobj, gene = gene)  
    }
  }
  if (scale.difference){
    mat.scale <- scalematrix(fit[gene, ,drop=F])
  } else {
    max <- apply(abs(fit[gene, ,drop=F]), 1, max)
    mat.scale <- fit[gene, ,drop=F]/max
  }
  
  if (method == 'kmeans'){
    set.seed(12345)
    # 
    if (k.auto){
      clu <- mykmeans(mat.scale, maxclunum = 20)$cluster
    } else {
      clu <- kmeans(mat.scale, k, iter.max = 1000)$cluster
    }
  } else if (method == 'hierarchical') {
    clu <- cutree(hclust(dist(mat.scale)), k = k)
  }
  
  # order clusters by genes' earliest max expr position
  
  v <- sapply(unique(clu), function(i){
    ap <- which(colMeans(mat.scale[names(clu)[clu==i], -ncol(mat.scale), drop=FALSE]) * colMeans(mat.scale[names(clu)[clu==i], -1, drop = FALSE]) < 0)
    ap[which.min(abs(ap-ncol(mat.scale)/2))]
  })
  names(v) <- unique(clu)
  
  corv <- apply(mat.scale,1,cor,1:ncol(mat.scale))
  corv <- tapply(corv,list(clu),mean)
  corv <- corv[names(v)]
  # self study
  v[corv < 0] <- ncol(mat.scale)-v[corv < 0]
  v <- v * (2*(corv > 0)-1)
  
  trans <- cbind(as.numeric(names(sort(v))),1:length(v))
  n <- names(clu)
  clu <- trans[match(clu,trans[,1]),2]
  
  names(clu) <- n
  
  clu2 <- paste0(clu, ';',rowMeans(fit[names(clu), , drop=F]) > 0)
  uclu2 <- sort(unique(clu2))
  clu2 <- match(clu2,uclu2)
  names(clu2) <- n
  return(clu2)  
}

# testobj = Res
# gene = diffgene
# type = 'variable'
# k=3
# k.auto = FALSE
# k=5
# method = 'kmeans'
# scale.difference = F
# 
# gene = names(DDGType)[!DDGType %in% c('nonDDG', 'meanSig')]
# type = 'variable'
# k=k
# scale.difference = scale.difference
# method = method
# k.auto = k.auto


clusterGene <- function(testobj, gene, type = 'variable', k.auto = FALSE,  k=5, method = 'kmeans', scale.difference = F){
  ## k.auto: if FALSE (default), use input k value. If TRUE, automatically select number of clusters. 
  if (toupper(type) == 'TIME'){
    clu <- cluster_gene(testobj = testobj, 
                        gene = gene,
                        k = k,
                        type = type,
                        method = method)
  } else if (toupper(type) == 'VARIABLE'){
    if ('DDGType' %in% names(testobj)){
      DDGType <- testobj$DDGType
    } else {
      DDGType <- getDDGType(testobj)
    }
    
    clu <- cluster_gene(testobj, gene = names(DDGType)[!DDGType %in% c('nonDDG', 'meanSig')], type = 'variable', k=k, scale.difference = scale.difference, method = method, k.auto = k.auto)
    design = testobj$design
    cellanno = testobj$cellanno
    meandiff <- sapply(c(0,1), function(i){
      as <- rownames(design[design[,2]==i, ])
      if ('expr' %in% names(testobj)){
        rowMeans(testobj$expr[names(DDGType)[DDGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1], drop = FALSE])
      } else {
        rowMeans(testobj$expr.ori[names(DDGType)[DDGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1], drop = FALSE])
      }
      
    }, simplify = FALSE)
    meandiff = do.call(cbind, meandiff)
    large0 <- rownames(meandiff)[meandiff[,1] >= meandiff[,2]]
    large1 <- rownames(meandiff)[meandiff[,1] < meandiff[,2]]
    
    clu2 <- rep(max(clu)+1, length(large1))
    names(clu2) <- large1
    clu3 <- rep(max(clu)+2, length(large0))
    names(clu3) <- large0
    clu = c(clu, clu2, clu3)  
  }
  clu <- clu[gene]
  return(clu)
}

