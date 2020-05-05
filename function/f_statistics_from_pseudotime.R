f_statistics_from_pseudotime <-  function(order.matrix, group.v1, group.v2){
  order = order.matrix
  colnames(order) <- c('Pseudotime', 'Cell', 'Patient')
  g1 = group.v1
  g2 = group.v2
  rm(order.matrix)
  tmin <- as.numeric(min(order$Pseudotime))
  tmax <- as.numeric(max(order$Pseudotime))
  res <- t(sapply(c(g1,g2), function(i){
    density(order[order$Patient == i, "Pseudotime"], from = tmin, to = tmax)$y
  }))
  rownames(res) <- c(g1,g2)
  colnames(res) <- density(order[order$Patient == g1[1], "Pseudotime"], from = tmin, to = tmax)$x
  ### spline fit
  trainX =   density(order[order$Patient == g1[1], "Pseudotime"], from = tmin, to = tmax)$x
  para = get_spline_coefficient(res, trainX, fit.min = tmin, fit.max = tmax)
  f <- ANOVA_model_f_stat(para, c(rep(1, length(g1)), rep(0, length(g2))))
  return(f)
  
}

