signal <- as.numeric(commandArgs(trailingOnly = T)[[1]])
method <- as.character(commandArgs(trailingOnly = T)[[2]])
# signal = 1
# method = 'EM'
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca')
rdir <- './simu/testvar/addMultiSignalUsingExpr/result/'
ddir <- './data/simu/testvar/addMultiSignalUsingExpr/'
fn <- paste0(rdir, method,'/', signal, '.rds')
print(fn)
if (file.exists(fn)) break 
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS(paste0('./data/simu/testtime/poolSampleSignal/null/pseudotime.rds'))
selgene <- readRDS(paste0(ddir, 'selgene/selgene.rds'))
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)

### two group along pseudotime
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  ### prepare data [Note: in two group senario, these should be different: design, v, cellWeights]
  expr <- readRDS(paste0(ddir, 'count/', signal, '.rds'))
  expr <- expr[rowMeans(expr > 0) > 0.01, ]
  design = matrix(c(1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  pdt <- data.frame(curve1 = pseudotime[,2], curve2 = pseudotime[,2])
  rownames(pdt) <- pseudotime[,1]
  pdt = pdt[colnames(expr), ]
  v <- (cellanno$sample %in% paste0('BM',c(1,2,5,6)) + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(counts)
  ### run test
  set.seed(12345)
  sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE)
  saveRDS(sce, paste0(rdir, method,'/', signal, '_sce.rds'))
  Final <- list()
  for (TestType in (c('diffEndTest', 'patternTest', 'earlyDETest'))){
    print(TestType)
    if (grepl('diffEndTest', TestType)){
      Res <- diffEndTest(sce)  
    } else if (grepl('patternTest', TestType)){
      Res <- patternTest(sce)  
    } else if (grepl('earlyDETest', TestType)){
      Res <- earlyDETest(sce, knots = c(1,2), global = TRUE, pairwise = TRUE)
    }
    res <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
    row.names(res) <- row.names(Res)
    res <- res[order(res[,3], -res[,1]), ]
    Final[[TestType]] <- list(res = res,
                              sensfdr = c(method, AreaUnderSensFdr(SensFdr(TruePositive = selgene, Statistics=res))))
  }
  saveRDS(Final, fn)  
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
  testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design,ncores=8, permuiter=100, type = 'Variable', demean = TRUE)
  saveRDS(testres, fn)  
  # ### calculate auc, fdr.diff
  # res <- data.frame(adj.P.Val = testres$fdr, foldchange = testres$foldchange, stringsAsFactors = F)
  # rownames(res) <- names(testres$fdr)
  # res <- res[order(res[,1], -res[,2]),,drop = FALSE]
  # testres[['sensfdr']] <-  c(method, AreaUnderSensFdr(SensFdr(TruePositive = selgene, Statistics = res)))
  # saveRDS(testres, paste0(rdir, method,'/', signal, '_perf.rds'))  
}

if (method == 'EM_NOT_centered'){
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
  testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design,ncores=4, permuiter=100, type = 'Variable', demean = FALSE)
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


