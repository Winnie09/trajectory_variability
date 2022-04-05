getCovariateGroupDiff <- function(testobj,
                                  gene, 
                                  reverse = FALSE) {
  ## testobj: output object from testpt().
  ## gene: a character vector of genes.
  ## output: a gene by pseudotime matrix. Entries are group difference w.r.t the variable, i.e., the unit-covariate incremental difference.
  ## reverse: logitcal. FALSE (default) when group 1 - group 0.  TRUE when group 0 - group 1. 
  ## (deprecated) variable: A character denoting the covariate for population fit. It should be one of the colnames in the design matrix.
  knotnum = testobj$knotnum[gene]
  pseudotime = seq(1, max(testobj$pseudotime))
  beta <- lapply(gene, function(g) {
    if (reverse){
      - testobj$parameter[[g]]$beta[c(seq(1, knotnum[g]+4), seq((testvar-1)*(knotnum[g] + 4)+1, testvar*(knotnum[g] + 4)))] ###
    } else {
      testobj$parameter[[g]]$beta[c(seq(1, knotnum[g]+4), seq((testvar-1)*(knotnum[g] + 4)+1, testvar*(knotnum[g] + 4)))] ###
    }
  })
  names(beta) <- gene
  
  philist <-
    lapply(min(knotnum):max(knotnum), function(num.knot) {
      if (num.knot == 0) {
        # phi <- cbind(1, bs(pseudotime))
        phi <- bs(pseudotime, intercept = TRUE)
      } else {
        knots = seq(min(pseudotime), max(pseudotime), length.out = num.knot + 2)[2:(num.knot +
                                                                                      1)]
        # phi <- cbind(1, bs(pseudotime, knots = knots))
        phi <- bs(pseudotime,knots = knots, intercept = TRUE)
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


