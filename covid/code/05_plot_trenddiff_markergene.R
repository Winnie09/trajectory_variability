
Res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/trenddiff.rds')
names(Res)
g <- c('ID1:ENSG00000125968', 'CCR10:ENSG00000184451', 'CXCL8:ENSG00000169429', 'CD9:ENSG00000010278', 'TMEM140:ENSG00000146859')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/'
pdf(paste0(pdir, 'trenddiff_markergene_sample_curve.pdf'), width = 14, height = 8)
plotGene(Res, g)
dev.off()
