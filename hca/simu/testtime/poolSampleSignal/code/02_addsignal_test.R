# --------------
# global setting
# --------------
clusterType <- as.numeric(commandArgs(trailingOnly = T)[[1]])
pctGene <- as.numeric(commandArgs(trailingOnly = T)[[2]])
method <- as.character(commandArgs(trailingOnly = T)[[3]])
# clusterType <- 1
# pctGene <- 1
# method <- 'tradeSeq'

print(method)
print(paste0('clusterType', clusterType, '_', pctGene))

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
ddir <- './hca/data/simu/testtime/poolSampleSignal/'
rdir <- './hca/simu/testtime/result/addsignal/'
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))

# ------------
# prepare data
# ------------
pt <- readRDS('./hca/data/simu/testtime/poolSampleSignal/null/pseudotime.rds')
pseudotime <- pt[,2]
names(pseudotime) <- pt[,1]

# ----------------------------------
# test  [one group along pseudotime]
# ----------------------------------
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  expr <- readRDS(paste0(ddir, 'count/clusterType', clusterType, '_', pctGene, '.rds'))
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'intercept')
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
  saveRDS(sce, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_sce.rds'))
  Final <- list()
  for (TestType in (c('startVsEndTest', 'associationTest'))){
    print(TestType)
    if (grepl('startVsEndTest', TestType)){
      Res <- startVsEndTest(sce)
    } else if (grepl('associationTest', TestType)){
      Res <- associationTest(sce)
    } 
    res <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
    row.names(res) <- row.names(Res)
    res <- res[order(res[,3], -res[,1]), ]
    Final[[TestType]] <- res
  }
  saveRDS(Final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}

if (method == 'EM_centered'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.1, ]
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  # design = cbind(1,design)
  testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time')
  saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}


if (method == 'EM_NOT_centered'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.1, ]
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  # design = cbind(1,design)
  testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time', demean = FALSE)
  saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}


if (method == 'tscan'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  expr <- expr[, names(pseudotime)]
  res <- TSCAN_time(expr=expr,pseudotime=pseudotime)
  saveRDS(res, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}

if (method == 'monocle2'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  expr <- expr[, names(pseudotime)]
  res <- monocle2_time(expr=expr,pseudotime=pseudotime)
  saveRDS(res, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}

if (method == 'monocle3'){
  library(spdep)
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  expr <- expr[, names(pseudotime)]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/poolSampleSignal/null/hsc_mep_ery_integrated_pca.rds')
  res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
  saveRDS(res, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}


