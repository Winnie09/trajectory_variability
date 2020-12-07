# --------------
# global setting
# --------------
dataType <- as.numeric(commandArgs(trailingOnly = T)[[1]])
method <- as.character(commandArgs(trailingOnly = T)[[2]])
# dataType <- 1
# method <- 'EM_NOT_centered'
print(paste0('dataType', dataType))
print(method)

suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))
ddir <- here('hca/data/simu/testtime/addMultiSignalUsingExpr/')
rdir <- here('hca/simu/testtime/result/addsignal/')
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)

# ------------
# prepare data
# ------------
pt <- readRDS(here('hca/data/simu/testtime/poolSampleSignal/null/pseudotime.rds'))
pseudotime <- pt[,2]
names(pseudotime) <- pt[,1]

# ----------------------------------
# test  [one group along pseudotime]
# ----------------------------------
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  expr <- readRDS(paste0(ddir, 'count/', dataType, '.rds'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  
  ##
  pdt <- data.frame(curve1 = pseudotime, curve2 = pseudotime)
  rownames(pdt) <- names(pseudotime)
  pdt = pdt[colnames(expr), ]
  
  v <- (cellanno$sample %in% paste0('BM',seq(1,8)) + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(expr)
  
  set.seed(12345)
  sce <- fitGAM(counts = round(expr), pseudotime = pdt, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE,parallel=F)
  # sce <- fitGAM(counts = round(expr), pseudotime = pdt, cellWeights = cellWeights,
  #               nknots = 6, verbose = FALSE,parallel=TRUE, BPPARAM = MulticoreParam(2))
  saveRDS(sce, paste0(rdir, method,'/', dataType,'_sce.rds'))
  res <- list()
  for (TestType in (c('startVsEndTest', 'associationTest'))){
    print(TestType)
    if (grepl('startVsEndTest', TestType)){
      Res <- startVsEndTest(sce)
    } else if (grepl('associationTest', TestType)){
      Res <- associationTest(sce)
    } 
    resdf <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
    row.names(resdf) <- row.names(Res)
    resdf <- resdf[order(resdf[,3], -resdf[,1]), ]
    res[[TestType]] <- resdf
  }
  
}

if (method == 'EM_centered'){
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.1, ]
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  res <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time', demean = TRUE, return.all.data = TRUE)
}


if (method == 'EM_NOT_centered'){
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.1, ]
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  res <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time', demean = FALSE, return.all.data = TRUE)
}


if (method == 'tscan'){
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- expr[, names(pseudotime)]
  res <- TSCAN_time(expr=expr,pseudotime=pseudotime)
}

if (method == 'monocle2'){
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- expr[, names(pseudotime)]
  res <- monocle2_time(expr=expr,pseudotime=pseudotime)
}

if (method == 'monocle3'){
  library(spdep)
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- expr[, names(pseudotime)]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/poolSampleSignal/null/hsc_mep_ery_integrated_pca.rds')
  res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
}

saveRDS(res, paste0(rdir, method,'/', dataType,'.rds'))  




