Res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/trenddiff.rds')
names(Res)

m <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/PRJCA002413_pbmc/data/proc/pt/expr.rds')

Res[['expr.demean']] = Res$expression
Res[['expr.ori']] = m[rownames(Res$expression), colnames(Res$expression)]

g <- c('ID1:ENSG00000125968', 'CCR10:ENSG00000184451', 'CXCL8:ENSG00000169429', 'CD9:ENSG00000010278', 'TMEM140:ENSG00000146859')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/'
pdf(paste0(pdir, 'trenddiff_markergene_sample_curve.pdf'), width = 14, height = 8)
plotGene(Res, g)
dev.off()

g = names(sort(Res$fdr[Res$fdr < 0.05]))
pdf(paste0(pdir, 'trenddiff_markergene_original_expr.pdf'), width = 25, height = 20)
plotGene(Res, g, original.expr = T, plot.point = T, point.alpha = 0.05, point.size = 0.05)
dev.off()
