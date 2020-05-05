#source('/Users/wenpinhou/Dropbox/trajectory_variability/function/get_spline_coefficient.R')
get_t_stat <- function(order, g1, g2){
  ## calculate the t-statistics of two groups in cell composition along a pseudotime trajectory.
  ## output is a vector of 14 t-statistics of 14 spline coefficients.
  ## order: a dataframe, 1st col is Pseudotime, 2nd col is Cell name. Cell name is in a format of "sampleName:cellName".
  ## g1: a vector of sample names in group 1.
  ## g2: a vector of sample names in group 2
  allp <- gsub(':.*','', order$Cell)  
  tmin = min(order$Pseudotime)
  tmax = max(order$Pseudotime)
  ## kernel density
  res <- t(sapply(c(g1,g2), function(i){
    density(order[allp == i, "Pseudotime"], from = tmin, to = tmax)$y
  }))
  rownames(res) <- c(g1,g2)
  ## spline fitting
  trainX =   density(order[allp == g1[1], "Pseudotime"], from = tmin, to = tmax)$x
  para = get_spline_coefficient(res, trainX, fit.min = tmin, fit.max = tmax)
  colnames(para) = paste0('coef',seq(1,ncol(para)))
  ## t statistics
  m1 = t(para[g1,])
  m2 = t(para[g2,])
  t <- sapply(rownames(m1), function(gene){
          v1 = m1[gene,]
          v2 = m2[gene,]
          (mean(v1) - mean(v2))/(sqrt(var(v1)/length(v1) + var(v2)/length(v2)))
  })
  t 
}

