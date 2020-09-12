get_population_fit <- function(testobj, variable = 'condition', value = NULL, gene){
  design = testobj$design
  pseudotime = testobj$pseudotime
  knotnum = testobj$knotnum
  pseudotime = pseudotime[order(pseudotime)]
  
  beta <- lapply(testobj$parameter[gene], function(i){
    i$beta
  })
  names(beta) <- gene
  if (sum(design[,1])!=nrow(design)) print("The first column of design matrix should be all 1s (intercept)!")
  colnames(design)[1] <- 'intercept'
  design <- unique(design[,c('intercept', variable)])

  rownames(design) <- paste0(variable, '_', unique(design[, variable]))

  x <- sapply(row.names(design),function(i) {
        kronecker(diag(knotnum[[gene]] + 4), design[i,])
      },simplify = F)
  maxknotallowed <- 30
  philist <- lapply(0:maxknotallowed,function(num.knot) {
    if (num.knot==0) {
      phi <- cbind(1,bs(pseudotime))
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      phi <- cbind(1,bs(pseudotime,knots = knots))
    }
  })
  names(philist) <- as.character(0:maxknotallowed)

  phi <- philist[[as.character(knotnum[gene])]]
  beta <- beta[[gene]]

  fit <- lapply(x, function(i){
    phi %*% t(i) %*% beta
  })
  names(fit) <- names(x)
  return(fit)
}
