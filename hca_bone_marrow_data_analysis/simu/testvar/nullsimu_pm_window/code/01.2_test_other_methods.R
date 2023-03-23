library(data.table)
library(limma)
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
### prepare data
expr <- readRDS(paste0(ddir, 'saverlog_pm.rds'))
expr <- expr[rowMeans(expr>0)>0.01, ]
design = readRDS(paste0(ddir, 'design.rds'))
cellanno = readRDS(paste0(ddir, 'cellanno_pm.rds'))
pt = readRDS(paste0(ddir, 'pseudotime_pm.rds'))
### run test

method = as.character(commandArgs(trailingOnly = T)[[1]][1])
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/', method)
dir.create(rdir, showWarnings = F, recursive = T)
fn = paste0(rdir, '/res.rds')
print(fn)


if (method == 'tscan'){
  branch = sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res = TSCAN_group(expr, pt, branch)
  saveRDS(res, fn)
}

if (method == 'meandiff'){
  agg <- vapply(unique(cellanno[,2]), function(i)
    rowMeans(expr[, cellanno[,2] == i, drop = FALSE], na.rm = TRUE), 
    numeric(nrow(expr)))
  agg <- agg[, paste0('BM', 1:8)]
  res <- topTable(eBayes(lmFit(agg, design)),n=nrow(agg),coef=2)
  saveRDS(res, fn)
}
