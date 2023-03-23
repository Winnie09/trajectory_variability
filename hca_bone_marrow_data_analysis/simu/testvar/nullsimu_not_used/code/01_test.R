suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
setwd(here())
# setwd(here('hca','simu','testvar','addMultiSignalUsingExpr','data'))
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
rdir <- 'hca/simu/testvar/nullsimu/result/'
ddir <- 'hca/simu/testvar/nullsimu/data/'

method <- as.character(commandArgs(trailingOnly = T)[[1]][1])
print(method)
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)
fn <- paste0(rdir, method, '/res.rds')
print(fn)

### load saver, count matrix, and pseudotime
saver = readRDS(paste0(ddir, 'saverlog_pm.rds'))
cnt = readRDS(paste0(ddir, 'cnt_pm.rds'))
cellanno = readRDS(paste0(ddir, 'cellanno_pm.rds'))
design = readRDS(paste0(ddir, 'design.rds'))
pt = readRDS(paste0(ddir, 'pseudotime_pm.rds'))



if (method == 'Lamian'){
  res <- testpt(expr=saver, cellanno=cellanno, pseudotime=pt, design=design, ncores=2, permuiter=100, test.type = 'Variable', demean = FALSE)
  saveRDS(res, fn)
}

if (method == 'tscan'){
  branch = sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res = TSCAN_group(saver, pt, branch)
  saveRDS(res, fn)
}

if (method == 'meandiff'){
  agg <- vapply(unique(cellanno[,2]), function(i)
    rowMeans(saver[, cellanno[,2] == i, drop = FALSE], na.rm = TRUE), 
    numeric(nrow(saver)))
  agg <- agg[, paste0('BM', 1:8)]
  res <- topTable(eBayes(lmFit(agg, design)),n=nrow(agg),coef=2)
  saveRDS(res, fn)
}


