library(tradeSeq)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/cnt_pm.rds')
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/'
pt = readRDS(paste0(ddir, 'pseudotime_pm.rds'))
d <- d[,names(pt)]
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/design.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/cellanno_pm.rds')
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]

x <- design[match(cellanno[,2],rownames(design)),2]
x <- factor(x)
sce <- fitGAM(counts = d, pseudotime = matrix(1:ncol(d),ncol=1), cellWeights = matrix(1,nrow=ncol(d),ncol=1),condition=x, verbose = TRUE)
saveRDS(sce,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/condiments/sce.rds')
cond_genes <- conditionTest(sce)
saveRDS(cond_genes,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/condiments/cond_genes.rds')

