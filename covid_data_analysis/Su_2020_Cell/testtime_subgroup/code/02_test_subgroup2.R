method <- as.character(commandArgs(trailingOnly = T)[[1]])
library(here)
setwd(here())
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_4levels.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]

source('function/01_function.R')
rdir <- paste0('covid/Su_2020_Cell/testtime_subgroup/result/', method, '/')
dir.create(rdir,recursive = T)
g1 <- unlist(sapply(unique(design[,2]),function(i) {
  tmp <- rownames(design[design[,2]==i,])
  set.seed(12345)
  sample(tmp,length(tmp)/2)
}))
g2 <- setdiff(rownames(design),g1)

##### select group
cellid <- cellanno[cellanno[,2] %in% g2,1]
expr=expr[,cellid]
cellanno=cellanno[cellid,]
pseudotime=pt[cellid]
design=design[g2,1,drop=F]
savefn <- paste0(rdir, 'testtime_res_subgroup2.rds')

# ----------------------------------
# test  [one group along pseudotime]
# ----------------------------------
if (method == 'EM_pm'){
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, test.type='Time', ncores = 16, demean = FALSE, test.method = 'permutation', ncores.fit = 8)  
}

if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  cnt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/data/count.rds')
  cnt <- cnt[rownames(expr), cellid]
  ##
  pdt <- data.frame(curve1 = pseudotime, curve2 = pseudotime)
  rownames(pdt) <- names(pseudotime)
  pdt = pdt[colnames(expr), ]
  
  v <- rep(1, nrow(cellanno))
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(expr)
  
  set.seed(12345)
  sce <- fitGAM(counts = round(expr), pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=F)
  saveRDS(sce, paste0(rdir, 'testtime_res_subgroup2_sce.rds'))

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

if (method == 'EM_chisq'){
  res <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design = design, test.type = 'Time',ncores = 48, permuiter=1000, test.method = 'chisq')
}

if (method == 'tscan'){
  res <- TSCAN_time(expr=expr,pseudotime=pseudotime)
}

if (method == 'monocle2'){
  res <- monocle2_time(expr=expr,pseudotime=pseudotime)
}

if (method == 'monocle3'){
  library(spdep)
  int <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/data/CD8integrate.rds')
  pca <- int@reductions$pca@cell.embeddings
  pca <- pca[cellid, ]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/poolSampleSignal/null/hsc_mep_ery_integrated_pca.rds')
  res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
}
saveRDS(res, savefn)


