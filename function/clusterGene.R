clusterGene <- function(testobj, 
                        gene,
                        k,
                        type = 'time',
                        method = 'kmeans'){
  ## testobj: output object from testpt().
  ## gene: a character vector of genes.
  ## k: number of clusters.
  ## type: "time" or "variable". A character denoting the population fit is on "time" or "variable".  Default is "time".
  
  if (type == 'time'){
    if ('populationFit' %in% names(testobj)) {
      fit <- testobj$populationFit
    } else {
      fit <- getPopulationFit(Res, gene = gene, type = type)
    }
  } else if (type == 'variable'){
    if ('covariateGroupDiff' %in% names(testobj)){
      fit <- testobj$covariateGroupDiff
    } else{
      fit <- getCovariateGroupDiff(testobj = testobj, gene = gene)  
    }
  }
  mat.scale <- scalematrix(fit[gene, ,drop=F])
  
  if (method == 'kmeans'){
      set.seed(12345)
      clu <- kmeans(mat.scale, k, iter.max = 1000)$cluster
  }
  ## order clusters by genes' max expr position
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



