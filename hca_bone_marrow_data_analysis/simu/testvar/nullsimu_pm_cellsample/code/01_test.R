suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
setwd(here())
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
rdir <- 'hca/simu/testvar/nullsimu/result/'

method <- as.character(commandArgs(trailingOnly = T)[[1]])

print(method)
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)
fn <- paste0(rdir, method, '/res.rds')
print(fn)

### load saver, count matrix, and pseudotime
pt <- readRDS('hca/simu/testvar/nullsimu/data/pseudotime_pm.rds')
saverlog <- readRDS('hca/simu/testvar/nullsimu/data/saverlog_pm.rds')
cellanno <- readRDS('hca/simu/testvar/nullsimu/data/cellanno_pm.rds')
design <- readRDS('hca/simu/testvar/nullsimu/data/design.rds')

expr <- saverlog[, names(pt)]
expr <- expr[rowMeans(expr>0.01)>0.01, ]


if (method == 'EM_pm'){
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, ncores=34, permuiter=100, test.type = 'Variable', demean = FALSE, test.method = 'permutation')
  saveRDS(res, fn)
}

if (method == 'EM_chisq'){
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, ncores=34, permuiter=100, test.type = 'Variable', demean = FALSE, test.method = 'chisq')
  saveRDS(res, fn)
}

if (method == 'tscan'){
  branch = sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res = TSCAN_group(expr, pt, branch, parallel = T)
  saveRDS(res, fn)
}

if (method == 'limma'){
  agg <- vapply(unique(cellanno[,2]), function(i)
    rowMeans(expr[, cellanno[,2] == i, drop = FALSE], na.rm = TRUE), 
    numeric(nrow(expr)))
  agg <- agg[, paste0('BM', 1:8)]
  res <- topTable(eBayes(lmFit(agg, design)),n=nrow(agg),coef=2)
  saveRDS(res, fn)
}

