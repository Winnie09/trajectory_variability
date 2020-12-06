get_population_fit <- function(testobj,
                               variable = 'condition',
                               value = NULL,
                               gene) {
  library(splines)
  ## if variable = NA, then return population fit (a vector) for constant test (test on time)
  ## if variable = a character (one of the colname in design matrix), then return population fit for all levels of that character.
  design = testobj$design
  pseudotime = testobj$pseudotime
  knotnum = testobj$knotnum
  pseudotime = pseudotime[order(pseudotime)]
  
  beta <- lapply(testobj$parameter[gene], function(i) {
    i$beta
  })
  names(beta) <- gene
  if (sum(design[, 1]) != nrow(design))
    print("The first column of design matrix should be all 1s (intercept)!")
  colnames(design)[1] <- 'intercept'
  if (is.na(variable)) {
    design = design[, 1, drop = FALSE]
  } else {
    design <- unique(design[, c('intercept', variable)])
    rownames(design) <-
      paste0(variable, '_', unique(design[, variable]))
  }
  
  
  x <- sapply(row.names(design), function(i) {
    kronecker(diag(knotnum[gene] + 4), design[i, , drop = FALSE]) ###
  }, simplify = FALSE)
  maxknotallowed <- 30
  philist <- lapply(0:maxknotallowed, function(num.knot) {
    if (num.knot == 0) {
      phi <- cbind(1, bs(pseudotime))
    } else {
      knots = seq(min(pseudotime), max(pseudotime), length.out = num.knot + 2)[2:(num.knot +
                                                                                    1)]
      phi <- cbind(1, bs(pseudotime, knots = knots))
    }
  })
  names(philist) <- as.character(0:maxknotallowed)
  
  phi <- philist[[as.vector(as.character(knotnum[gene]))]]
  beta <- beta[[gene]]
  
  if (is.na(variable)) {
    fit <- t(phi %*% t(x[[1]]) %*% beta)[1,]
  } else {
    fit <- lapply(x, function(i) {
      phi %*% t(i) %*% beta
    })
    names(fit) <- names(x)
  }
  return(fit)
}

