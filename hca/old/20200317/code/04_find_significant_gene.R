setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA')
## find variable genes
# mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/proc/matrix/saver.rds')
# order = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/order.rds')
# mat = mat[, order$Cell]
# ap = sub('_.*','',colnames(mat))
# gs <- sapply(unique(ap), function(s){
#   smat =  mat[,grepl(s,colnames(mat))]
#   smat = smat[rowMeans(smat>0.1)>0.1,]
#   cv <- sapply(rownames(smat), function(i){
#     sd(smat[i,])/mean(smat[i,])
#   })
#   print(quantile(cv,0.5,na.rm=T))
#   print(median(cv,na.rm=T))
#   names(cv[cv>quantile(cv,0.5,na.rm=T)])
# })
# 
# g = unique(unlist(gs))
# saveRDS(g, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/variable_genes_cv0.5.rds')
# g = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/variable_genes_0.1.rds')
g = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/variable_genes_cv0.5.rds')

#### read files
allfn = list.files('./ery')
lenv = NULL
for (fn in allfn){
  lenv[fn] = length(list.files(paste0('./ery/', fn)))
}
lenv = lenv[lenv==9]

pval = 0
for (fn in names(lenv)){
  p = readRDS(paste0('./ery/', fn,'/pvalue.rds'))
  pval = p + pval
}
# pval = pval/(length(lenv[as.numeric(names(lenv))<10])*100 + length(lenv[as.numeric(names(lenv))>10])*20)
pval = pval/(100*length(lenv))

int = intersect(g,rownames(pval))
pval = pval[int,]

for (i in 1:ncol(pval)){
  pval[, i] = p.adjust(pval[, i], method='fdr')
}
sig <- (pval<0.05) + 0
genes <- rownames(sig)[which(rowSums(sig) > 0)]
df = data.frame(gene = rownames(sig),num.significant.base = rowSums(sig,na.rm = T), stringsAsFactors = F)
df = df[order(df[,2],decreasing = T),]
saveRDS(df,'./ery/significant_genes.rds')


