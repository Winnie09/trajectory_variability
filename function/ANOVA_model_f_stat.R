ANOVA_model_f_stat <- function(mat, group){
  ## this function compare two model fitting to see if the groups parameters difference is zero
  ## input mat: sample by num.parameters
  ## input group: a numeric vector indicating two groups, length(group)==nrow(mat). For example, group = c(1,1,0,0) indicates four samples and the first two belong to first group and last two belong to the other group.
  ## output: a f statistics calculating from full and null model RRS
  ## Full model
  group = as.numeric(as.vector(group))
  sse1 <- sum(sapply(seq(1,ncol(mat)), function(i){
    v = mat[,i]
    x = matrix(c(rep(1,nrow(mat)), group), ncol=2)
    betatmp = solve(t(x)%*%x) %*% t(x)%*% v
    sum((v - x %*% betatmp)^2)
  }))
  
  # Null model
  sse2 <- sum(sapply(seq(1,ncol(mat)), function(i){
    v = mat[,i]
    x = matrix(rep(1,nrow(mat)), ncol=1)
    betatmp = solve(t(x)%*%x) %*% t(x)%*% v
    sum((v - x %*% betatmp)^2)
  }))
  ((sse2-sse1)/14)/(sse1/(length(mat)-28))
}

