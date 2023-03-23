# --------------
# global setting
# --------------
method <- as.character(commandArgs(trailingOnly = T)[[1]])
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
setwd(here())
ddir <- 'hca/real/testtime_twogroup/data/monocyte/'
rdir <- paste0('hca/real/testtime_twogroup/result/', method, '/monocyte/3478/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

# ------------
# prepare data
# ------------
expr = readRDS(paste0(ddir, 'expr2.rds'))
cellanno = readRDS(paste0(ddir, 'cellanno2.rds'))
rownames(cellanno) <- cellanno[,1]
pseudotime = readRDS(paste0(ddir, 'pseudotime2.rds'))
design = readRDS(paste0(ddir, 'design2.rds'))
expr <- expr[rowMeans(expr > 0.01) > 0.01, ]
cnt <- readRDS('hca/data/proc/matrix/count.rds')
cnt <- cnt[, colnames(expr)]
cnt <- cnt[rowMeans(cnt > 0.01) > 0.01, ]
# ----------------------------------
# test  [one group along pseudotime]
# ----------------------------------
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  expr <- cnt
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
  saveRDS(sce, paste0(rdir, 'sce.rds'))
  
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

if (method == 'EM_pm'){
  system.time({
    res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=24, test.type = 'Time', demean = FALSE, test.method = 'permutation', ncores.fit = 48)
  })
}

if (method == 'EM_chisq'){
  system.time({
    res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=24, test.type = 'Time', demean = FALSE, test.method = 'chisq', ncores.fit = 48)
  }) 
}

if (method == 'tscan'){
  res <- TSCAN_time(expr=expr,pseudotime=pseudotime)
}

if (method == 'monocle2'){
  res <- monocle2_time(expr=expr,pseudotime=pseudotime)
}

if (method == 'monocle3'){
  library(spdep)
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- expr[, names(pseudotime)]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/poolSampleSignal/null/hsc_mep_ery_integrated_pca.rds')
  res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
}
saveRDS(res, paste0(rdir, 'testtime_res.rds'))


