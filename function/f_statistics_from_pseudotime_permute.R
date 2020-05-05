f_statistics_from_pseudotime_permute <- function(order.matrix, group.v1, group.v2, num.permute=1000){
  f_permute <- sapply(seq(1, num.permute), function(myseed){
    g = c(group.v1, group.v2)
    set.seed(myseed)
    id = sample(1:nrow(order.matrix),nrow(order.matrix), replace=T)  ## cell 
    set.seed(myseed)
    g = sample(g) ## group
    g1 = g[seq(1,length(group.v1))]
    g2 = g[seq((length(group.v1)+1), length(g))]
    f <- f_statistics_from_pseudotime(order.matrix[id,], g1, g2)
  })
  return(f_permute)
}

