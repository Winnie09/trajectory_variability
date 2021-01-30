f_statistics_from_gene_permute <- function(expr,order.matrix, group.v1, group.v2, num.permute=3){
  expr <- expr[,order.matrix[,2]]
  library(parallel)
  f_permute <- mclapply(seq(1, num.permute), function(myseed){
    print(myseed)
    g = c(group.v1, group.v2)
    set.seed(myseed)
    id = sample(1:nrow(order.matrix),nrow(order.matrix), replace=T)  ## cell 
    set.seed(myseed)
    g = sample(g) ## group
    g1 = g[seq(1,length(group.v1))]
    g2 = g[seq((length(group.v1)+1), length(g))]
    f <- f_statistics_from_gene(expr[,id],order.matrix[id,,drop=F], g1, g2)
  },mc.cores=20)
  f_permute = do.call(cbind, f_permute)
  colnames(f_permute) = paste0('pm',seq(1,num.permute))
  return(f_permute)
}

