ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
### prepare data
expr <- readRDS(paste0(ddir, 'saverlog_pm.rds'))
expr <- expr[rowMeans(expr>0)>0.01, ]
design = readRDS(paste0(ddir, 'design.rds'))
cellanno = readRDS(paste0(ddir, 'cellanno_pm.rds'))
pt = readRDS(paste0(ddir, 'pseudotime_pm.rds'))
### run test
system.time({
  testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, ncores=16, test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'permutation')
})
saveRDS(testres, paste0(rdir, 'testvar_res.rds')) 


