getPopulationFit <- function(testobj,
                             gene = NULL,
                             type = 'time',
                             num.timepoint = 1e3){
  library(splines)
  ## if type = 'time', then return population fit (a vector for a gene; or a gene by num.cell matrix) for constant test (test on time)
  ## if type = 'variable', then return population fit for all levels of that character (a matrix, columns are population fit for each level in the variabel). 
  ## gene: a vector of gene names. 
  design = testobj$design[, c(1, testobj$testvar)] ## design for multi
  knotnum = testobj$knotnum
  pseudotime = testobj$pseudotime
  pseudotime = pseudotime[order(pseudotime)]
  pt <- round(seq(1, max(pseudotime), length.out = min(num.timepoint, max(pseudotime)))) ## downsample
  testvar = testobj$testvar
  type <- toupper(type)
  if (sum(design[, 1]) != nrow(design)){
    print("The first column of design matrix should be all 1s (intercept)! Using the first column as the variable column ...")
    design = cbind(intercept = 1, design)
  }
  colnames(design)[1] <- 'intercept'
  if (is.null(gene)) gene <- rownames(testobj$statistics)
  if (type == 'TIME') {
    design = design[, 1, drop = FALSE]
  } else {
    variable = colnames(design)[2]
    design <- unique(design[, c(colnames(design)[1], variable)])
    rownames(design) <- paste0(variable, '_', unique(design[, variable]))
  }
  
  fitlist <- lapply(gene, function(g){
    print(g)
    # beta <- lapply(testobj$parameter[g], function(i) {
    #   i$beta
    # })
    # names(beta) <- g
    tmp = matrix(testobj$parameter[[g]]$beta, ncol = knotnum[g]+4)
    beta = as.vector(tmp[c(1,testvar), ]) ### subset the beta values of the intercept and the test covariate for multi
    x <- sapply(row.names(design), function(i) {
      kronecker(diag(knotnum[g] + 4), design[i, , drop = FALSE]) ###
    }, simplify = FALSE)
    
    
    if (knotnum[g] == 0) {
      # phi <- cbind(1, bs(pt))
      phi <- bs(pt, intercept = TRUE)
    } else {
      knots = seq(min(pt), max(pt), length.out = knotnum[g] + 2)[2:(knotnum[g] + 1)]
      # phi <- cbind(1, bs(pt, knots = knots))
      phi <- bs(pt,knots = knots, intercept = TRUE)
    }
    
    if (exists('variable')) {
      fit <- lapply(x, function(i) {
        if (ncol(phi) == nrow(i)){
          phi %*% i %*% beta
        } else {
          phi %*% t(i) %*% beta
        }
      })
    } else {
      i = x[[1]]
      if (ncol(phi) == nrow(i)){
        fit <- phi %*% i %*% beta
      } else {
        fit <- phi %*% t(i) %*% beta
      }
    }
    return(fit)
  })
  
  names(fitlist) <- gene
  if (type == 'VARIABLE'){
    fitres <- lapply(names(fitlist[[1]]), function(i){
      tmp <- t(sapply(fitlist, function(j){
        j[[i]]
      }))
    })  
    names(fitres) <- names(fitlist[[1]])
  } else if (type == 'TIME'){
    fitres <- t(do.call(cbind, fitlist))
    rownames(fitres) <- gene
    if (ncol(Res$expr) == ncol(fitres)) colnames(fitres) <- colnames(Res$expr)
  }
  return(fitres)
}


