setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA')
## find variable genes
mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/proc/matrix/saver.rds')
order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/order.rds')
mat = mat[, order$Cell]
ap = sub('_.*','',colnames(mat))
gs <- sapply(unique(ap), function(s){
  smat =  mat[,grepl(s,colnames(mat))]
  smat = smat[rowMeans(smat>0.1)>0.1,]
  cv <- sapply(rownames(smat), function(i){
    sd(smat[i,])/mean(smat[i,])
  })
  print(quantile(cv,0.75,na.rm=T))
  print(median(cv,na.rm=T))
  names(cv[cv>quantile(cv,0.75,na.rm=T)])
})
  
g = unique(unlist(gs))
saveRDS(g, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/variable_genes_0.1.rds')

