library(here)
setwd(here())
path = as.character(commandArgs(trailingOnly = T)[[1]])
ddir <- paste0('hca/real/build_from_tree_variability/result/', path, '/')
rdir <- paste0('hca/real/testvar/result/tradeSeq/', path, '/gender/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

source('./function/01_function.R')
expr = readRDS(paste0(ddir, 'input_expr.rds'))
cnt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/matrix/count.rds')
expr = cnt[intersect(rownames(expr), rownames(cnt)), colnames(expr)]
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)

### two group along pseudotime
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(slingshot))
suppressMessages(library(tradeSeq))
### prepare data [Note: in two group senario, these should be different: design, v, cellWeights]

pdt <- data.frame(curve1 = pseudotime, curve2 = pseudotime)
rownames(pdt) <- names(pseudotime)
pdt = pdt[colnames(expr), ]

v <- (cellanno$sample %in% rownames(design)[design[,2]==0] + 0)
v <- ifelse(v==1, 0.99, 0.01)

cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
rownames(cellWeights) <- colnames(expr)
### run test
set.seed(12345)
sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE,parallel=TRUE,BPPARAM=BiocParallel::MulticoreParam(workers = 15))
saveRDS(sce, paste0(rdir, 'sce.rds'))

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
  Final[[TestType]] <-  res
}
saveRDS(Final, paste0(rdir, 'testvar_res.rds'))  

