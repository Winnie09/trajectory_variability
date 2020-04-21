f_statistics_from_gene <- function(expression,order.matrix, group.v1, group.v2){
  ## expression: gene by cell matrix
  ## order.matrix: a dataframe, columns are (ordered) Pseudotime, Cell, Patient
  ## group.v1: a character vector indicating sample names in group 1
  ## group.v2: a character vecgtor indicating sample names in group 2
  order = order.matrix
  colnames(order) <- c('Pseudotime', 'Cell', 'Patient')
  g1 = group.v1
  g2 = group.v2
  rm(order.matrix)
  tmin <- as.numeric(min(order$Pseudotime))
  tmax <- as.numeric(max(order$Pseudotime))
  allg <- c(g1,g2)
  #expression <- expression[,order[,2]]
  res <- lapply(allg,function(sg) get_spline_coefficient(expression[,order[,3]==sg],order[order[,3]==sg,1],tmin,tmax))
  names(res) <- allg
  sapply(row.names(res[[1]]),function(g) {
    para <- t(sapply(res,function(i) i[g,]))
    row.names(para) <- allg
    ANOVA_model_f_stat(para, c(rep(1, length(g1)), rep(0, length(g2))))
  })
}

