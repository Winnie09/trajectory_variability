setwd('/home/whou10/scratch16/whou10/trajectory_variability/')
d <- readRDS('hca/simu/testvar/nullsimu_pm_window/data/data/cnt_pm.rds')
ddir <- 'hca/simu/testvar/nullsimu_pm_window/data/data/'
psn = readRDS(paste0(ddir, 'pseudotime_pm.rds'))
expr <- d[,names(psn)]
design <- readRDS('hca/simu/testvar/nullsimu_pm_window/data/data/design.rds')
cellanno <- readRDS('hca/simu/testvar/nullsimu_pm_window/data/data/cellanno_pm.rds')
cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]

suppressMessages(library(SingleCellExperiment))
suppressMessages(library(slingshot))
suppressMessages(library(tradeSeq))
### prepare data [Note: in two group senario, these should be different: design, v, cellWeights]

pdt <- data.frame(curve1 = psn, curve2 = psn)
rownames(pdt) <- names(psn)
pdt = pdt[colnames(expr), ]

v <- (cellanno$sample %in% rownames(design)[design[,2]==0] + 0)
v <- ifelse(v==1, 0.99, 0.01)

cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
rownames(cellWeights) <- colnames(expr)
### run test
set.seed(12345)
sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE,parallel=FALSE)
rdir <- 'hca/simu/testvar/nullsimu_pm_window/result/'
dir.create(paste0(rdir, 'tradeseq'), recursive = TRUE)
saveRDS(sce, paste0(rdir, 'tradeseq/sce.rds'))

Final <- list()
for (TestType in c('diffEndTest', 'patternTest', 'earlyDETest')){
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
  Final[[TestType]] <- res
}
saveRDS(Final, paste0(rdir, 'tradeseq/res.rds'))  


