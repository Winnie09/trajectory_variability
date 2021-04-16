cluster_gene <- function(testobj, 
                         gene,
                         k,
                         type = 'time',
                         method = 'kmeans', 
                         scale.difference = F){
  ## testobj: output object from testpt().
  ## gene: a character vector of genes.
  ## k: number of clusters. 2 are for meanSig, and the remaining k-2 clusters are for DDG of NOT meanSig.
  ## type: "time" or "variable". A character denoting the population fit is on "time" or "variable".  Default is "time". ## method: 'kmeans'(default) or 'hierarchical'.
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
    clu <- kmeans(mat.scale, k, iter.max = 1000)$cluster
  } else if (method == 'hierarchical') {
    clu <- cutree(hclust(dist(mat.scale)), k = k)
  }
  # order clusters by genes' max expr position
  v <- sapply(unique(clu), function(i) {
    tmp <- fit[clu == i,]
    mean(apply(tmp, 1, which.max))
  })
  names(v) <- unique(clu)
  trans <- cbind(as.numeric(names(v)),rank(v))
  n <- names(clu)
  clu <- trans[match(clu,trans[,1]),2]
  names(clu) <- n
  return(clu)  
}


clusterGene <- function(testobj, gene, type = 'variable', k=5, method = 'kmeans', scale.difference = F){
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
    clu <- cluster_gene(testobj, gene = names(DDGType)[!DDGType %in% c('nonDDG', 'meanSig')], type = 'variable', k=k-2, scale.difference = scale.difference, method = method)
    
    
    design = testobj$design
    cellanno = testobj$cellanno
    meandiff <- sapply(c(0,1), function(i){
      as <- rownames(design[design[,2]==i, ])
      if ('expr' %in% names(testobj)){
        rowMeans(testobj$expr[names(DDGType)[DDGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1]])
      } else {
        rowMeans(testobj$expr.ori[names(DDGType)[DDGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1]])
      }
        
    })
    large0 <- rownames(meandiff)[meandiff[,1] >= meandiff[,2]]
    large1 <- rownames(meandiff)[meandiff[,1] < meandiff[,2]]
    
    clu2 <- rep(k-1, length(large0))
    names(clu2) <- large0
    clu3 <- rep(k, length(large1))
    names(clu3) <- large1
    clu = c(clu, clu2, clu3)  
  }
  clu <- clu[gene]
  return(clu)
}

