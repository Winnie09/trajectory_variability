# --------------
# global setting
# --------------
method <- as.character(commandArgs(trailingOnly = T)[[1]])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './hca/simu/testtime/result/nullsimu/'
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')

# ------------
# prepare data
# ------------
# expression
expr <- readRDS('./hca/data/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/proc/ct/sc.rds')
ct <- ct[colnames(expr)]
expr <- expr[, ct %in% c('HSC','MEP','Ery')]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
# cellanno
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
# design
design <- cbind(rep(1,8), c(1,1,0,0,1,1,0,0))
dimnames(design) = list(paste0('BM',1:8), c('intercept', 'group'))
# pseudotime [reorder cells for each sample]
pt <- readRDS('./hca/data/simu/testtime/null/pseudotime.rds')
pt <- data.frame(pt, sample = sub(':.*','', pt$cell), stringsAsFactors = FALSE)
ptlist <- lapply(unique(pt[,3]), function(p){
  tmp <- pt[as.character(pt[,3])==p, ]
  set.seed(12345)
  tmp[, 'pseudotime'] <- sample(tmp[, 'pseudotime'])
  tmp
})
pt <- do.call(rbind, ptlist)
pseudotime <- pt[,2]
names(pseudotime) <- pt[,1]

# ----
# test
# ----
if (method == 'EM_SelectKnots'){
  testres <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=8, permuiter=100, type = 'Time')
  saveRDS(testres, paste0(rdir, method, '/testres.rds'))
}
  
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  ##  one grouop
  counts <- round(exp(expr + 1))
  pdt <- data.frame(curve1 = pseudotime, curve2 = pseudotime)
  rownames(pdt) <- names(pseudotime)
  pdt = pdt[colnames(expr), ]
  
  v <- (cellanno$sample %in% paste0('BM',seq(1,8)) + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(counts)
  
  set.seed(12345)
  sce <- fitGAM(counts = counts, pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE)
  saveRDS(sce, paste0(rdir, method,'/sce.rds'))
  
  #######
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
    saveRDS(Final, paste0(rdir, method,'/testres.rds'))  
}

if (method == 'tscan'){
  expr <- expr[, names(pseudotime)]
  res <- TSCAN_time(expr=expr,pseudotime=pseudotime)
  saveRDS(res, paste0(rdir, method,'/testres.rds'))  
}

if (method == 'monocle2'){
  expr <- expr[, names(pseudotime)]
  res <- monocle2_time(expr=expr,pseudotime=pseudotime)
  saveRDS(res, paste0(rdir, method,'/testres.rds'))
}

if (method == 'monocle3'){
  expr <- expr[, names(pseudotime)]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/null/hsc_mep_ery_integrated_pca.rds')
  res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
  saveRDS(res, paste0(rdir, method,'/testres.rds'))
}


