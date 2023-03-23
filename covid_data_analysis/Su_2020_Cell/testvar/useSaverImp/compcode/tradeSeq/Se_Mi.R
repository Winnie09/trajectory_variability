setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
method <- 'tradeSeq'

source('./function/01_function.R')
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Se_Mi.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]
#cnt <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/su/count.rds')
#cnt <- cnt[rownames(expr),colnames(expr)]

ser <- readRDS('covid/Su_2020_Cell/data/CD8integrate.rds')
cnt <- as.matrix(ser@assays$RNA@counts)
cnt <- cnt[intersect(rownames(cnt),rownames(expr)), colnames(expr)]

rdir <- paste0('/home-4/zji4@jhu.edu/scratch/diffpt/su/compres/tradeSeq/Se_Mi/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)


### two group along pseudotime
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(slingshot))
suppressMessages(library(tradeSeq))
### prepare data [Note: in two group senario, these should be different: design, v, cellWeights]

pdt <- data.frame(curve1 = pt, curve2 = pt)
rownames(pdt) <- names(pt)
pdt = pdt[colnames(expr), ]

v <- (cellanno$sample %in% rownames(design)[design[,2]==0] + 0)
v <- ifelse(v==1, 0.99, 0.01)

cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
rownames(cellWeights) <- colnames(expr)
### run test
set.seed(12345)
sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE,parallel=TRUE)
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
  Final[[TestType]] <- res
}
saveRDS(Final, paste0(rdir, 'testvar_res.rds'))  
