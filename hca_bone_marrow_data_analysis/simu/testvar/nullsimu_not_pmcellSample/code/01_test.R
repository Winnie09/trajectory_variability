suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
setwd(here())
# setwd(here('hca','simu','testvar','addMultiSignalUsingExpr','data'))
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
rdir <- 'hca/simu/testvar/nullsimu/result/'

method <- as.character(commandArgs(trailingOnly = T)[[1]])
print(method)
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)
fn <- paste0(rdir, method, '/res.rds')
print(fn)

### load saver, count matrix, and pseudotime
saverlog <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','hsc_mep_ery_saver.rds'))
pt <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','pseudotime.rds'))
expr <- saverlog[, pt[,1]]
expr <- expr[rowMeans(expr>0)>0.01, ]
design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
pseudotime <- pt[, 2]
names(pseudotime) <- pt[, 1]


if (method == 'EM_pm'){
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=16, permuiter=100, test.type = 'Variable', demean = FALSE, test.method = 'permutation')
  saveRDS(res, fn)
}

if (method == 'EM_chisq'){
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=24, test.type = 'Variable', demean = FALSE, test.method = 'chisq')
  saveRDS(res, fn)
}


if (method == 'tscan'){
  branch = sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res = TSCAN_group(expr, pseudotime, branch)
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




