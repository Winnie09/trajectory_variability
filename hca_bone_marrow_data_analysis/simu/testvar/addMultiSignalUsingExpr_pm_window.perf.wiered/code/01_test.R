signal <- as.character(commandArgs(trailingOnly = T)[[1]])
method <- as.character(commandArgs(trailingOnly = T)[[2]])
# signal = 1
# method = 'EM'
print(signal)
print(method)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca')
rdir <- 'simu/testvar/addMultiSignalUsingExpr_pm_window/result/'
ddir <- 'simu/testvar/addMultiSignalUsingExpr_pm_window/data/data/'
fn <- paste0(rdir, method,'/', signal, '.rds')
print(fn)
if (file.exists(fn)) break 
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS('simu/testvar/nullsimu_pm_window/data/data/pseudotime_pm.rds')
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)

### two group along pseudotime
if (method == 'Lamian.chisq'){
  print(method)
  ### prepare data
  expr <- readRDS(paste0(ddir, 'saver/', signal, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr>0)>0.01, ]
  
  design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  ### run test
  testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design,ncores=10, test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'chisq')
  saveRDS(testres, fn)  
}

if (method == 'Lamian.pm'){
  print(method)
  ### prepare data
  expr <- readRDS(paste0(ddir, 'saver/', signal, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr>0)>0.01, ]
  
  design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  ### run test
  system.time({
    testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=4, test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'permutation')
  })
  saveRDS(testres, fn)  
}


warnings()

# expr: saver imputed matrix
# count: count matrix
# pseudotime: numeric vector (1,2,3,4....) with names same as colnames(expr)
# branch: 0,1 vector indicating whether each cell is from group 1 or 2, can get from as.numeric(sub(':.*','',colnames(expr)) %in% paste0('BM',c(1,2,5,6)))
# cell_coords: the pca you sent me, only use the first 4 (if correct) dimensions




