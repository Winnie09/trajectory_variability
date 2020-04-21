get_spline_fit <- function(trainData, trainX, fit.min, fit.max, num.base=10, remove.correlated=TRUE){
  ## trainData: gene by cell matrix
  ## trainX: vector, e.g. pseudotime of cells, same length as ncol(trainData)
  ## fit.min, fit.max: scalar, fitting range
  ## output: gene by spline fitting matrix
  knots = seq(fit.min,fit.max,length.out=num.base+2)[2:(num.base+1)]
  library(splines)
  base = cbind(1,bs(trainX,knots = knots))
  if (remove.correlated){
    colidx = NULL
    for (ii in seq(2,ncol(base))){
      if (length(unique(base[,ii])) == 1) colidx = c(colidx, ii)
    }
    if (length(colidx)) base = base[,-colidx]
  }
  coef = t(chol2inv(chol(crossprod(base))) %*% t(base) %*% t(trainData) )
  pred = t(base %*% t(coef))
}

