regression_ANOVA_f_stat <- function(mat, group){
  ## this function perform linear regression on two groups of samples and then ANOVA test
  ## input: mat, numSample by numCoefficient matrix
  ## input: group, a numeric vector denoting groups. length(group)==nrow(mat). e.g. group <- c(1,1,1,0,0,0) denotes six samples and the first three ones are in group 1 while the other in group 2
  ## output: ANOVA test f statis
  group = as.numeric(as.vector(group))
  beta <- sapply(seq(1,ncol(mat)), function(i){
  v = mat[,i]
  x = matrix(c(rep(1,nrow(mat)), group), ncol=2)
  betatmp = solve(t(x)%*%x) %*% t(x)%*% v
  })
  
  v1 = beta[1,]
  v2 = beta[2,]
  f <- get_ANOVA_f_stat(beta[1,], beta[2,]+beta[1,])
  return(f)
}

