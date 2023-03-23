rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/'
Res <- readRDS(paste0(rdir, 'temra_trenddiff.rds'))
names(Res)
res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange, pvalue = Res$pvalue)
res <- res[res[,1] < 0.05, ]
rownames(res) <- sub(':.*', '', rownames(res))
write.csv(res, paste0(rdir, 'temra_trenddiff.csv'))

dg = rownames(res[res[,1]<0.05, ])
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/'

pdf(paste0(pdir, 'temra_trenddiff_gene_sample_curve.pdf'), width = 14, height = 8)
plotGene(Res, dg)
dev.off()

pdf(paste0(pdir, 'temra_trenddiff_gene_sample_curve_type_demean.pdf'), width = 14, height = 8)
plotGene(Res, dg, variable= 'type', variable.text = c('Healthy', 'COVID-19'))
dev.off()

pdf(paste0(pdir, 'temra_trenddiff_gene_sample_curve_type_ori.pdf'), width = 14, height = 8)
plotGene(Res, dg, original.expr = TRUE, variable= 'type', variable.text = c('Healthy', 'COVID-19'))
dev.off()




