rm(list=ls())
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
order = readRDS('./hca/result/ery/order.rds')
mat = readRDS('./hca/data/HCA/proc/matrix/saver.rds')
mat = mat[,order$Cell]
mat <- mat[rowMeans(mat>0.01)>0.1, ]
source('./function/01_function.R')
order = data.frame(order, Patient = gsub('_.*','', order$Cell))
ap = as.character(unique(order$Patient))
g1 = ap[grepl('female', ap)]
g2 = ap[grepl(':male', ap)]
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
# source('/Users/wenpinhou/Dropbox/resource/function.R')
eg <- sapply(ap, function(p){
  print(p)
  tmat <- mat[,grepl(p, colnames(mat))]
  rownames(tmat)[rowMeans(tmat>0.01)>0.1]
})
eg <- unique(unlist(eg))
mat = mat[eg,]
vg <- sapply(ap, function(p){
  print(p)
  tmat <- mat[,grepl(p, colnames(mat))]
  vg <- findVariableGene(tmat, num.gene = NULL ,plot.statistics=TRUE, plot.dir = paste0('./hca/geneexpr/plot.mac/',p,'/'))
})
vg <- unique(unlist(vg))
mat = mat[vg,]
b = readRDS('./hca/geneexpr/result/f_statistics_from_gene_gender.rds')
a = readRDS('./hca/geneexpr/result/f_statistics_from_lowExprGene_gender_permute_new1e4.rds')
a = a[rownames(mat),]
b = b[rownames(a)]
pval <- sapply(seq(1,nrow(a)), function(i){
  sum(a[i,]>b[i])/ncol(a)
})
names(pval) = rownames(a)
fdr = p.adjust(pval,method='fdr')

# ag <- names(sort(pval)[1:16])
ag <- names(sort(b, decreasing=TRUE)[1:16])
library(ggplot2)
library(gridExtra)
plist <- list()
for (g in ag){
  print(g)
  pd1 = mat[g, grepl('female', colnames(mat))]
  pd1 = data.frame(Expr=pd1, Cell=names(pd1), Patient = gsub('_.*','',names(pd1) ), Gender='Female')
  pd2 = mat[g, grepl(':male', colnames(mat))]
  pd2 = data.frame(Expr=pd2, Cell=names(pd2), Patient = gsub('_.*','',names(pd2) ), Gender='Male')
  pd = rbind(pd1, pd2)
  pd = cbind(pd, Pseudotime = order[match(pd$Cell, order$Cell),'Pseudotime'])
  linedlist <- lapply(unique(pd$Patient), function(p){
    tmat = mat[g,grepl(p,colnames(mat)),drop=F]
    trainX = order$Pseudotime[grepl(p,colnames(mat))]
    pred <- get_spline_fit(tmat, trainX=seq(1,ncol(tmat)), fit.min=min(order$Pseudotime), fit.max=max(order$Pseudotime))
    tmpdf <- data.frame(Expr=pred[1,], Pseudotime=trainX, Patient=p, Gender=ifelse(grepl('female',p),'female','male'))
  })
  ld = do.call(rbind, linedlist)
  plist[[g]] <- ggplot() + geom_point(data=pd, aes(x=Pseudotime, y=Expr, color=Patient), alpha=.1, size=.2)  + 
    geom_line(data=ld, aes(x=Pseudotime, y=Expr, color=Patient),alpha=1, size=.5) +
    theme_classic() + ggtitle(paste0(sub(':.*','',g),',p=', round(pval[g],3),',f=',round(b[g],2))) + theme(legend.position = 'none') + scale_color_manual(values=c(rep('darkblue',4),rep('orange',4)))
  }
pdf('./hca/geneexpr/plot.mac/gender_diff_gene_top_f.pdf',width=12,height=9)
# pdf('./hca/geneexpr/plot.mac/gender_diff_gene_top_pval.pdf',width=12,height=9)
grid.arrange(grobs=plist,nrow=4)
dev.off()


############# plot order permutation result
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
allg = sub(':.*','',names(sort(pval)))
str(allg)

v1 <- sapply(seq(1,length(allg)), function(i){
  sum(allg[seq(1,i)] %in% u1)
})
  
v2 <- sapply(seq(1,length(allg)), function(i){
  sum(allg[seq(1,i)] %in% u2)
})
v1_pm <- sapply(seq(1,1e2), function(myseed){
  set.seed(myseed)
  w1 = sample(allg, length(u1))
  v1 <- sapply(seq(1,length(allg)), function(i){
    sum(allg[seq(1,i)] %in% w1)
  })
})
rownames(v1_pm) <- paste0('top',seq(1,nrow(v1_pm)))
saveRDS(v1_pm, './hca/geneexpr/result/geneset_same_length_as_chrX_gene_pm_mean_order.rds')
v1_pm <- rowMeans(v1_pm)
v2_pm <- sapply(seq(1, 1e2), function(myseed){
  print(myseed)
  set.seed(myseed)
  w2 = sample(allg, length(u2))  
  v2 <- sapply(seq(1,length(allg)), function(i){
    sum(allg[seq(1,i)] %in% w2)
  })
})  
rownames(v2_pm) <- paste0('top',seq(1,nrow(v2_pm)))
saveRDS(v2_pm, './hca/geneexpr/result/geneset_same_length_as_chrY_gene_pm_mean_order.rds')
v2_pm <- rowMeans(v2_pm)
df = data.frame(chrX=v1, chrY=v2, chrX_pm = v1_pm, chrY_pm = v2_pm, order = seq(1,length(v1)))
saveRDS(df, './hca/geneexpr/result/df_chrX_chrY_pm_order.rds')
mat <- NULL
for (i in 1:4) {
  mat <- rbind(mat,data.frame(v=df[,i],order=df[,5],type=colnames(df)[i]))
}
library(ggplot2)
pdf('./hca/geneexpr/plot.mac/chrX_chrY_order_compare_to_permutation.pdf', width=4, height=4)
ggplot(mat,aes(x=order,y=v,col=type, fill=type), alpha=.2) + geom_line() + xlim(c(0,30)) + ylim(c(0,10))+theme_classic()+ylab('number of ChrX/Y genes') + xlab('top n genes (ordered by increasing pvalue)')
dev.off()

## all chrX + chrY
u = unique(c(u1,u2))
v <- sapply(seq(1,length(allg)), function(i){
  sum(allg[seq(1,i)] %in% u)
})
v_pm <- sapply(seq(1, 1e2), function(myseed){
  print(myseed)
  set.seed(myseed)
  w = sample(allg, length(u))  
  v <- sapply(seq(1,length(allg)), function(i){
    sum(allg[seq(1,i)] %in% w)
  })
})  
rownames(v_pm) <- paste0('top',seq(1,nrow(v_pm)))
saveRDS(v_pm, './hca/geneexpr/result/geneset_same_length_as_chrXY_gene_pm_mean_order.rds')
v_pm = rowMeans(v_pm)

df = data.frame(chrXY = v, chrXY_pm = v_pm)
saveRDS(df, './hca/geneexpr/result/df_chrXY_pm_order.rds')
mat <- NULL
for (i in 1:2) {
  mat <- rbind(mat,data.frame(v=df[,i],order=seq(1,nrow(df)),type=colnames(df)[i]))
}

pdf('./hca/geneexpr/plot.mac/chrXY_order_compare_to_permutation.pdf', width=4, height=4)
ggplot(mat,aes(x=order,y=v,col=type, fill=type), alpha=.2) + geom_line() + xlim(c(0,30)) + ylim(c(0,10))+theme_classic()+ylab('number of ChrX/Y genes') + xlab('top n genes (ordered by increasing pvalue)')
dev.off()




