# --------------
# global setting
# --------------
library(here)
## setwd(here())
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
## setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
dataType <- as.numeric(commandArgs(trailingOnly = T)[[1]])
method <- as.character(commandArgs(trailingOnly = T)[[2]])
# dataType <- 1
# method <- 'EM_NOT_centered'
print(paste0('dataType', dataType))
print(method)

suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))
ddir <- here('hca/simu/testtime/addMultiSignalUsingExpr/data/')
rdir <- here('hca/simu/testtime/addMultiSignalUsingExpr/result/addsignal/')
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)

# ------------
# prepare data
# ------------
pt <- readRDS(here('hca/simu/testtime/poolSampleSignal/data/null/pseudotime.rds'))
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


if (method == 'EM_pm'){
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.1, ]
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  res <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design = design, test.type = 'Time',ncores = 48, permuiter=100, test.method = 'permutation')
}


if (method == 'EM_chisq'){
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.1, ]
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept'))
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  res <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design = design, test.type = 'Time',ncores = 48, permuiter=1000, test.method = 'chisq', cutoff = 1e-5, sd.adjust =2.6e-3)
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

if (method == 'pseudotimeDE'){ ## R4.0.2
  suppressMessages(library(SingleCellExperiment))
  suppressPackageStartupMessages(library(PseudotimeDE))
  suppressPackageStartupMessages(library(SingleCellExperiment))
  suppressPackageStartupMessages(library(slingshot))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(irlba))
  
  expr <- readRDS(paste0(ddir, 'count/', dataType, '.rds'))
  logexpr <- log2(expr + 1)
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  mysce <- SingleCellExperiment(list(counts=expr, logcounts = logexpr)) ## !!
  rd <- irlba::prcomp_irlba(t(logcounts(mysce)), scale. = FALSE)$x[, 1:2]
  
  reducedDims(mysce) <- SimpleList(PCA = rd)
  colData(mysce)$cl <- 1
  
  fit_ori <- slingshot(mysce, reducedDim = 'PCA', clusterLabels = "cl")
  ori_tbl <- tibble(cell = colnames(mysce), pseudotime = rescale(colData(fit_ori)$slingPseudotime_1))
  
  set.seed(123)
  ## Set the cores for parallelization. Note that mclapply doesnot work on Windows.
  options(mc.cores = 2)
  ## Number of subsmaples
  n = 100
  ## Ganerate random subsamples
  index <- mclapply(seq_len(n), function(x) {
    sample(x = c(1:dim(mysce)[2]), size = 0.8*dim(mysce)[2], replace = FALSE)
  })
  sub_tbl <- mclapply(index, function(x, sce) {
    sce <- sce[, x]
    rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
    reducedDims(sce) <- SimpleList(PCA = rd)
    
    fit <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "cl")
    tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
    
    ## Make sure the direction of pseudotime is the same as the original pseudotime
    merge.tbl <- left_join(tbl, ori_tbl, by = "cell")
    
    if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
      tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
    }
    tbl
  }, sce = mysce)
  
  system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = rownames(expr),
                                                   ori.tbl = ori_tbl,
                                                   sub.tbl = sub_tbl, ## To save time, use 100 subsamples
                                                   mat = mysce, ## You can also use a matrix or SeuratObj as the input
                                                   model = "nb",
                                                   mc.cores = 2))
}
saveRDS(res, paste0(rdir, method,'/', dataType,'.rds'))  





