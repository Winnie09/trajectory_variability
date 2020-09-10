Res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/trenddiff.rds')
names(Res)
res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange, pvalue = Res$pvalue)
res = res[order(res[,1], -res[,2]), ]
write.csv(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/trenddiff_full.csv')

Res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/meandiff.rds')
names(Res)
res = data.frame(fdr = Res$adj.P.Val, foldchange = Res$logFC, pvalue = Res$P.Value)
res = res[order(res[,1], -abs(res[,2])), ]
write.csv(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/meandiff_full.csv')

