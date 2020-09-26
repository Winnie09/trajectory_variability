Res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/trenddiff.rds')
names(Res)
res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange, pvalue = Res$pvalue)
res <- res[res[,1] < 0.05, ]
rownames(res) <- sub(':.*', '', rownames(res))
write.csv(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/trenddiff.csv')

dg = rownames(res[res[,1]<0.05, ])
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
expr = Res$expression
pred = Res$predict.values
rownames(expr) <- sub(':.*', '', rownames(expr))
Res$expression = expr
rownames(pred) = sub(':.*', '', rownames(pred))
Res$predict.values = pred

pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/'
pdf(paste0(pdir, 'trenddiff_gene_sample_curve.pdf'), width = 14, height = 8)
plotGene(Res, dg)
dev.off()


