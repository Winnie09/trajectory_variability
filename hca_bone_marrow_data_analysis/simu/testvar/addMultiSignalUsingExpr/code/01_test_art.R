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
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/null/pseudotime_pm_df.rds')
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)

### two group along pseudotime
if (method == 'EM_chisq'){
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
  expr <- expr[,names(pt)]
  
  #gid <- sample(setdiff(rownames(expr),selgene),5000)
  selgene <- readRDS('simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene.rds')
  gid <- setdiff(rownames(expr),selgene)
  selgene1 <- readRDS('simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds') #trend only
  selgene2 <- readRDS('simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds') #mean only
  selgene3 <- readRDS('simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds') #trend mean
#  gid <- c(gid,selgene2,selgene3)
  cid <- cellanno[cellanno[,2] %in% rownames(design)[design[,2]==1],1]
#  cid <- cid[1:length(cid) %in% sample(1:length(cid),length(cid)*0.5)]
  addmat <- matrix(rep(seq(1,length(cid))/length(cid)*3.3e-3,length(gid)),nrow=length(gid),byrow = T)
  expr[gid,cid] <- expr[gid,cid] + addmat
  #5e-3 too strong
#2e-3 too weak
#0.5 7e-3
#0.2 1.8e-2

  ### run test
  testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design,ncores=20, test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'chisq')
  saveRDS(testres, fn)  
}


if (method == 'EM_centered'){
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
  testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design,ncores=6, test.type = 'Variable', demean = TRUE)
  saveRDS(testres, fn)  
}

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
    testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, ncores=16, test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'permutation')
  })
  saveRDS(testres, fn)  
}

if (method == 'tscan'){
  print(method)
  expr <- readRDS(paste0(ddir, 'saver/', signal, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr>0)>0.01, ]
  design = matrix(c(1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  design = cbind(1,design)
  psn = pseudotime[,2]
  names(psn) = pseudotime[,1]
  branch = sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res = TSCAN_group(expr, psn, branch)
  saveRDS(res, fn)  
}

if (method == 'monocle2'){
  expr <- readRDS(paste0(ddir, 'saver/', signal, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr>0)>0.01, ]
  design = matrix(c(1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  design = cbind(1,design)
  psn = pseudotime[,2]
  names(psn) = pseudotime[,1]
  branch = sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res = monocle2_group(expr, psn, branch)
  saveRDS(res, fn)  
}

if (method == 'monocle3'){
  expr <- readRDS(paste0(ddir, 'saver/', signal, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr>0)>0.01, ]
  design = matrix(c(1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  
  design = cbind(1,design)
  psn = pseudotime[,2]
  names(psn) = pseudotime[,1]
  branch = sapply(1:nrow(cellanno), function(i) ((design[cellanno[i, 2], 2] == 1) + 0))
  res = monocle3_group(expr, branch)
  saveRDS(res, fn)  
}

warnings()

# expr: saver imputed matrix
# count: count matrix
# pseudotime: numeric vector (1,2,3,4....) with names same as colnames(expr)
# branch: 0,1 vector indicating whether each cell is from group 1 or 2, can get from as.numeric(sub(':.*','',colnames(expr)) %in% paste0('BM',c(1,2,5,6)))
# cell_coords: the pca you sent me, only use the first 4 (if correct) dimensions





