get_fitted_values <- function(testObj, Gene, Pseudotime=NULL){
  ### testObj: test result from testpt_Time
  library(splines)
  mat = matrix(0, nrow = length(Gene),ncol=length(Pseudotime),dimnames = list(Gene,NULL))
  num.knot = testObj$knotnum[Gene]
  if (is.null(Pseudotime)){
    Pseudotime = testObj$pseudotime
  }
  for (i in unique(num.knot)){
    x <- kronecker(diag(i + 4), 1)
    genes = names(num.knot[num.knot == i])
    if (i==0) {
      phi <- cbind(1,bs(Pseudotime))
    } else {
      knots = seq(min(Pseudotime),max(Pseudotime),length.out=i+2)[2:(i+1)]
      phi <- cbind(1,bs(Pseudotime,knots = knots))  
    }
    beta = sapply(genes,function(i) testObj$parameter[[i]]$beta)
    mat[genes,] <- t(phi %*% x %*% beta)
  } 
  return(mat)
}
