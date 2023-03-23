suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
setwd(here())
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
rdir <- 'hca/simu/testvar/nullsimu/result/'

method <- 'tradeSeq'

print(method)
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)
fn <- paste0(rdir, method, '/res.rds')
print(fn)

### load saver, count matrix, and pseudotime
psn <- readRDS('hca/simu/testvar/nullsimu/data/pseudotime_pm.rds')
expr <- readRDS('hca/simu/testvar/nullsimu/data/cnt_pm.rds')
cellanno <- readRDS('hca/simu/testvar/nullsimu/data/cellanno_pm.rds')
design <- readRDS('hca/simu/testvar/nullsimu/data/design.rds')

expr <- expr[, names(psn)]
expr <- expr[rowMeans(expr>0)>0.01, ]


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
              nknots = 6, verbose = FALSE,parallel=TRUE)

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
saveRDS(Final, fn)  
