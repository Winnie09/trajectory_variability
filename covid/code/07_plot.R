rdir <- '/dcl02/hongkai/data/whou/trajectory_variability/covid/result/'
Res <- readRDS(paste0(rdir, 'tex_trenddiff.rds'))
names(Res)
res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange, pvalue = Res$pvalue)
res <- res[res[,1] < 0.05, ]
rownames(res) <- sub(':.*', '', rownames(res))
write.csv(res, paste0(rdir, 'tex_trenddiff.csv'))

dg = rownames(res[res[,1]<0.05, ])
source('/dcl02/hongkai/data/whou/trajectory_variability/function/01_function.R')
pdir <- '/dcl02/hongkai/data/whou/trajectory_variability/covid/plot/'
pdf(paste0(pdir, 'tex_trenddiff_gene_sample_curve.pdf'), width = 14, height = 8)
plotGene(Res, dg)
dev.off()

