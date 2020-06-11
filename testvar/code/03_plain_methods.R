geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
# geneProp <- 0.05
# addSignalType <- 'linear'
# addSignalPara <-  1
method = as.character(commandArgs(trailingOnly = T)[[4]])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testvar/result/plain/'
datadir = './testvar/data/plain/'

fn = paste0(rdir, method, '/', addSignalType, '/', geneProp, '_', addSignalPara, '.rds')
if (file.exists(fn)) break

suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
source('./function/01_function.R')

d <- readRDS(paste0(datadir, addSignalType, '/', geneProp, '_', addSignalPara, '.rds'))
expr <- d$expr
selgene <- d$selgene
design <- d$design
dir.create(paste0(rdir, method, '/', addSignalType), recursive = T, showWarnings = F)
pseudotime <- readRDS('./testtime/data/data/null/pseudotime.rds')

cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)

if (method == 'EM_SelectKnots'){
  testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100)
  saveRDS(testres, paste0(sub('.rds','',fn), '_res.rds'))
  res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
  rownames(res) <- names(testres$fdr)
  res <- res[order(res[,1]),,drop=F]
  sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
  final <- list()
  final[['res']] <- res
  final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
  final[['perll']] <- testres$perll
  final[['knotnum']] <- testres$knotnum
  saveRDS(final, fn)  
}

if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  
  counts <- round(exp(expr + 1))
  pdt <- data.frame(curve1 = pseudotime[,2], curve2 = pseudotime[,2])
  rownames(pdt) <- pseudotime[,1]
  pdt = pdt[colnames(counts), ]
  
  v <- (cellanno$sample %in% paste0('BM',c(1,2,5,6)) + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(counts)
  
  set.seed(12345)
  sce <- fitGAM(counts = counts, pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE)
  saveRDS(sce, paste0(sub('.rds', '', fn), '_sce.rds'))
  
  Final <- list()
  for (TestType in (c('diffEndTest', 'patternTest', 'earlyDETest'))){
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
    sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
    final <- list()
    final[['res']] <- res
    final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
    Final[[TestType]] <- final
  }
  saveRDS(Final, fn)  
} 

