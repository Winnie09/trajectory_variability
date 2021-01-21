getPopulationFit <- function(testobj,
                             gene,
                             type = 'time'){
  library(splines)
  ## if type = 'time', then return population fit (a vector) for constant test (test on time)
  ## if type = 'variable', then return population fit for all levels of that character (a matrix, columns are population fit for each level in the variabel). 
  ## gene: a vector of gene names. 
  design = testobj$design
  pseudotime = testobj$pseudotime
  knotnum = testobj$knotnum
  pseudotime = pseudotime[order(pseudotime)]
  if (sum(design[, 1]) != nrow(design)){
    print("The first column of design matrix should be all 1s (intercept)! Using the first column as the variable column ...")
    design = cbind(intercept = 1, design)
  }
  colnames(design)[1] <- 'intercept'
  if (type == 'time') {
    design = design[, 1, drop = FALSE]
  } else {
    variable = colnames(design)[2]
     design <- unique(design[, c('intercept', variable)])
    rownames(design) <- paste0(variable, '_', unique(design[, variable]))
  }
  
  fitlist <- lapply(gene, function(g){
    beta <- lapply(testobj$parameter[g], function(i) {
      i$beta
    })
    names(beta) <- g
    x <- sapply(row.names(design), function(i) {
      kronecker(diag(knotnum[g] + 4), design[i, , drop = FALSE]) ###
    }, simplify = FALSE)
    
    if (knotnum[g] == 0) {
      phi <- cbind(1, bs(pseudotime))
    } else {
      knots = seq(min(pseudotime), max(pseudotime), length.out = knotnum[g] + 2)[2:(knotnum[g] + 1)]
      phi <- cbind(1, bs(pseudotime, knots = knots))
    }
    if (exists('variable')) {
      if (ncol(phi) == ncol(x[[1]])){
         fit <- t(phi %*% t(x[[1]]) %*% beta)[1,]
       } else {
         fit <- t(phi %*% x[[1]] %*% beta)[1,]
       }
    } else {
      fit <- lapply(x, function(i) {
        if (ncol(phi) == nrow(i)){
          phi %*% i %*% beta[[g]]
        } else {
          phi %*% t(i) %*% beta[[g]]
        }
      })
      names(fit) <- names(x)
    }
    return(fit)
  })
  names(fitlist) <- gene
  fitres <- lapply(names(fitlist[[1]]), function(i){
    tmp <- t(sapply(fitlist, function(j){
      j[[i]]
    }))
  })  
  names(fitres) <- names(fitlist[[1]])
  return(fitres)
}



