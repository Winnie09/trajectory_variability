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
# source('/Users/wenpinhou/Dropbox/resource/function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
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
b = readRDS('./hca/geneexpr/result/f_statistics_from_gene_age.rds')
a = readRDS('./hca/geneexpr/result/f_statistics_from_lowExprGene_age_permute_new1e4.rds')
a = a[rownames(mat),]
b = b[rownames(a)]
pval <- sapply(seq(1,nrow(a)), function(i){
  sum(a[i,]>b[i])/ncol(a)
})
names(pval) = rownames(a)
fdr = p.adjust(pval,method='fdr')
pdf('./hca/geneexpr/plot.mac/age_diff_f_p_fdr.pdf', width=7,height=4)
par(mfrow=c(1,2))
smoothScatter(pval~b, xlab='f statistics', ylab='p-value')
smoothScatter(fdr~b, xlab='f statistics', ylab='fdr')
dev.off()
ag <- names(sort(pval)[1:16])
# ag <- names(sort(b, decreasing=TRUE)[1:16])
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
# pdf('./hca/geneexpr/plot.mac/age_diff_gene_top_f.pdf',width=12,height=9)
pdf('./hca/geneexpr/plot.mac/age_diff_gene_top_pval.pdf',width=12,height=9)
grid.arrange(grobs=plist,nrow=4)
dev.off()



