signal <- as.character(commandArgs(trailingOnly = T)[[1]])
method <- as.character(commandArgs(trailingOnly = T)[[2]])
# signal = 1
# method = 'EM'
print(signal)
print(method)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca')
rdir <- 'simu/testvar/addMultiSignalUsingExpr/result/'
ddir <- 'simu/testvar/addMultiSignalUsingExpr/data/'
fn <- paste0(rdir, method,'/', signal, '.rds')
print(fn)
if (file.exists(fn)) break 
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/null/pseudotime_pm_df.rds')
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)

### two group along pseudotime

if (method == 'EM_pm'){
  print(method)
  ### prepare data
  expr <- readRDS(paste0(ddir, 'saver/', signal, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr>0)>0.01, ]
  
  design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  pt <- pseudotime[, 2]
  names(pt) <- pseudotime[, 1]
  ### run test
  system.time({
    testres <- testpt(expr=expr[1:10, ], cellanno=cellanno, pseudotime=pt, design=design, ncores=1, test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'permutation')
  })
  saveRDS(testres, fn)  
}

