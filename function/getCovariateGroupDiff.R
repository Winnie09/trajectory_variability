getCovariateGroupDiff <- function(testobj,
                               gene) {
  ## testobj: output object from testpt().
  ## gene: a character vector of genes.
  ## variable: A character denoting the covariate for population fit. It should be one of the colnames in the design matrix.
  ## output: a gene by pseudotime matrix. Entires are group difference w.r.t the variable, i.e., the unit-covariate incremental difference.
  knotnum = testobj$knotnum[gene]
  pseudotime = testobj$pseudotime
  beta <- lapply(testobj$parameter[gene], function(i) {
    i$beta
  })
  names(beta) <- gene
  philist <-
    lapply(min(knotnum):max(knotnum), function(num.knot) {
      if (num.knot == 0) {
        phi <- cbind(1, bs(pseudotime))
      } else {
        knots = seq(min(pseudotime), max(pseudotime), length.out = num.knot + 2)[2:(num.knot +
                                                                                      1)]
        phi <- cbind(1, bs(pseudotime, knots = knots))
      }
    })
  names(philist) <- as.character(min(knotnum):max(knotnum))
  
  fit <- lapply(gene, function(i) {
    id <- (1:(length(beta[[i]]) / 2)) * 2
    a <- philist[[as.character(knotnum[i])]] %*% beta[[i]][id]
  })
  fit <- do.call(cbind, fit)
  colnames(fit) <- gene
  return(t(fit))
}
