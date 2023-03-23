ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/lamian.chisq/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
### prepare data
expr <- readRDS(paste0(ddir, 'saverlog_pm.rds'))
expr <- expr[rowMeans(expr>0)>0.01, ]
design = readRDS(paste0(ddir, 'design.rds'))
cellanno = readRDS(paste0(ddir, 'cellanno_pm.rds'))
pt = readRDS(paste0(ddir, 'pseudotime_pm.rds'))
### run test

colnames(expr) <- paste0(sample(sub('_.*','',colnames(expr))),'_',sub('.*_','',colnames(expr)))
cellanno <- data.frame(cell=colnames(expr),sample=sub(':.*','',colnames(expr)),stringsAsFactors = F)
names(pt) <- colnames(expr)

expr[,cellanno[cellanno[,2] %in% rownames(design)[design[,2]==1],1]] <- expr[,cellanno[cellanno[,2] %in% rownames(design)[design[,2]==1],1]] + 2.9e-3
#2.8e-3 too small
#3e-3 too large

testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, ncores=20, test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'chisq')
sum(testres[[1]][,1] < 0.05)

saveRDS(testres, paste0(rdir, 'testvar_res.rds')) 

