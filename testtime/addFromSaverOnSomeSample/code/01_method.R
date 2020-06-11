geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
method <- as.character(commandArgs(trailingOnly = T)[[4]])
# geneProp <- 0.05
# addSignalType <- 'constant'
# addSignalPara <-  0.1
# method <- 'EM_SelectKnots'

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('./function/01_function.R')
datadir <- './testtime/data/data/'
rdir <- './testtime/result/'
if (!file.exists(paste0(rdir, method, '/', addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))){
  suppressMessages(library(parallel))
  suppressMessages(library(splines))
  suppressMessages(library(limma))
  suppressMessages(library(RColorBrewer))
  dir.create(paste0(rdir, method, '/', addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
  res <- readRDS(paste0(datadir,addSignalType,'/', geneProp,'_', addSignalPara,'.rds')) #####
  expr <- res$expr
  selgene <- res$selgene
  pseudotime <- res$pseudotime
  design <- res$design
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  
  ## two group along pseudotime
  if (grepl('tradeSeq', method)){
    suppressMessages(library(SingleCellExperiment))
    suppressMessages(library(slingshot))
    suppressMessages(library(tradeSeq))
    if (!file.exists(paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))){
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
      saveRDS(sce, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))
    } else {
      sce <- readRDS(paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))
    }
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
    saveRDS(Final, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  
  }
  
  if (method == 'EM_SelectKnots'){
    print(method)
    design = cbind(1,design)
    testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time')
    saveRDS(testres, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_testres.rds'))  
    res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
    rownames(res) <- names(testres$fdr)
    res <- res[order(res[,1]),,drop=F]
    sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
    final <- list()
    final[['res']] <- res
    final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
    final[['perll']] <- testres$perll
    final[['knotnum']] <- testres$knotnum
    saveRDS(final, paste0(rdir, method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  
  }
  rm(list=ls())
}


