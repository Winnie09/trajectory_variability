clusterType <- as.numeric(commandArgs(trailingOnly = T)[[1]])
pctGene <- as.numeric(commandArgs(trailingOnly = T)[[2]])
method <- as.character(commandArgs(trailingOnly = T)[[3]])
geneProp <- 0.05


setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('./function/01_function.R')
datadir <- './testtime/data/data/'
rdir <- './testtime/result/'
ddir <- './testtime/data/data/'
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)
pseudotime <- readRDS('./testtime/data/data/null/pseudotime.rds')
selgene = readRDS('./testtime/data/data/selgene/selgene.rds')
## one group along pseudotime
if (grepl('tradeSeq', method)){
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))
  expr <- readRDS(paste0(ddir, 'count/clusterType', clusterType, '_', pctGene, '.rds'))
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  
    ##
    pdt <- data.frame(curve1 = pseudotime[,2], curve2 = pseudotime[,2])
    rownames(pdt) <- pseudotime[,1]
    pdt = pdt[colnames(expr), ]
    
    v <- (cellanno$sample %in% paste0('BM',seq(1,8)) + 0)
    v <- ifelse(v==1, 0.99, 0.01)
    cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
    rownames(cellWeights) <- colnames(counts)
    
    set.seed(12345)
    sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
                  nknots = 6, verbose = FALSE,parallel=TRUE)
    saveRDS(sce, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_sce.rds'))
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
  saveRDS(Final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}

if (method == 'EM_SelectKnots'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  design = matrix(rep(1,8), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  expr <- log2(expr + 1)
  print(method)
  design = cbind(1,design)
  testres <- testpt(expr=expr,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100, type = 'Time')
  saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_testres.rds'))  
  res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
  rownames(res) <- names(testres$fdr)
  res <- res[order(res[,1]),,drop=F]
  # sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
  final <- list()
  final[['res']] <- res
  # final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
  final[['perll']] <- testres$perll
  final[['knotnum']] <- testres$knotnum
  saveRDS(final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}

if (method == 'tscan'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  psn <- as.numeric(pseudotime[,2])
  names(psn) <- pseudotime[,1]
  expr <- expr[, pseudotime[,2]]
  testres <- TSCAN_time(expr=expr,pseudotime=psn)
  saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_testres.rds'))  
  res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
  rownames(res) <- names(testres$fdr)
  res <- res[order(res[,1]),,drop=F]
  sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
  final <- list()
  final[['res']] <- testres
  final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
  saveRDS(final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}

if (method == 'monocle2'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  psn <- as.numeric(pseudotime[,2])
  names(psn) <- pseudotime[,1]
  expr <- expr[, pseudotime[,2]]
  testres <- monocle2_time(expr=expr,pseudotime=psn)
  saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_testres.rds'))  
  res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
  rownames(res) <- names(testres$fdr)
  res <- res[order(res[,1]),,drop=F]
  sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
  final <- list()
  final[['res']] <- testres
  final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
  saveRDS(final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}


if (method == 'monocle3'){
  expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
  psn <- as.numeric(pseudotime[,2])
  names(psn) <- pseudotime[,1]
  expr <- expr[, pseudotime[,2]]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/hsc_mep_ery_integrated_pca.rds')
  testres <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
  saveRDS(testres, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_testres.rds'))  
  res <- data.frame(adj.P.Val = testres$fdr, stringsAsFactors = F)
  rownames(res) <- names(testres$fdr)
  res <- res[order(res[,1]),,drop=F]
  sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
  final <- list()
  final[['res']] <- testres
  final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
  saveRDS(final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
}



