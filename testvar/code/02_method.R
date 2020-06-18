clusterType <- as.numeric(commandArgs(trailingOnly = T)[[1]])
pctGene <- as.numeric(commandArgs(trailingOnly = T)[[2]])
method <- as.character(commandArgs(trailingOnly = T)[[3]])

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testvar/result/'
ddir <- './testvar/data/data/'
fn <- paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds')
print(fn)

if (file.exists(fn)) break 

suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
source('./function/01_function.R')
pseudotime <- readRDS('./testtime/data/data/null/pseudotime.rds')
selgene <- readRDS('./testvar/data/data/selgene/selgene.rds')
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)


### two group along pseudotime
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  ### prepare data [Note: in two group senario, these should be different: design, v, cellWeights]
  expr <- readRDS(paste0(ddir, 'count/clusterType', clusterType, '_', pctGene, '.rds'))
  design = matrix(c(1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  pdt <- data.frame(curve1 = pseudotime[,2], curve2 = pseudotime[,2])
  rownames(pdt) <- pseudotime[,1]
  pdt = pdt[colnames(expr), ]
  v <- (cellanno$sample %in% paste0('BM',c(1,2,5,6)) + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- colnames(counts)
  ### run test
  set.seed(12345)
  sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE)
  saveRDS(sce, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_sce.rds'))
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

if (method == 'EM_SelectKnots'){
  ### prepare data
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  design = matrix(c(1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  expr <- log2(expr + 1)
  print(method)
  design = cbind(1,design)
  ### run test
  testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Variable')
  saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_testres.rds'))  
  ### calculate auc, fdr.diff
  res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
  rownames(res) <- names(testres$fdr)
  res <- res[order(res[,1]),,drop=F]
  sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
  testres[['sensfdr']] <-  c(method, AreaUnderSensFdr(sensfdr))
  saveRDS(testres, fn)  
}

warnings()
