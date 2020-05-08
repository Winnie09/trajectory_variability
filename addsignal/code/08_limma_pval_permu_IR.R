rm(list=ls())
geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
method <- as.character(commandArgs(trailingOnly = T)[[4]])
## be reminded that the pca is done on cv > 0.5 (for all samples) genes. You should redo it !!!
# geneProp <- 0.05
# addSignalType <- 'linear'
# addSignalPara <-  5
# method <- 'tradeSeq'

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
library(parallel)
library(splines)
library(limma)
suppressMessages(library(tradeSeq))
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)


expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
allp = sub(':.*','', colnames(expr))

dr <- readRDS('./addsignal/result/dr.rds')

library(TSCAN)
pseudotime <- TSCANorder(exprmclust(t(dr),reduce = F,clusternum=4,clustermethod='kmeans'),orderonly = T)
psn <- 1:length(pseudotime)
names(psn) <- pseudotime
order = data.frame(Pseudotime = seq(1, length(pseudotime)), Cell=pseudotime)


sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

## one group along pseudotime
rtgene <- sapply(c(1,0), function(i){
  print(i)
  selectsample <- rownames(design[design[,1]==i,,drop=F])
  print(str(selectsample))
  selectcell <- names(sample[sample %in% selectsample])
  res <- OneGroupDiffALongTime(expr = expr[, selectcell], sample=sample[selectcell], dr = dr[selectcell,], pseudotime=pseudotime[pseudotime %in% selectcell], num.base = 3) 
  rownames(res[res[,3]<0.05,])
})
expr <- expr[intersect(rtgene[[1]], rtgene[[2]]), ] ##  [1:2795, 1:13269]  ##  [1:6234, 1:13269]

set.seed(12345)
selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))

## two group along pseudotime
if (method == 'limma'){
  dir.create(paste0('./addsignal/result/limma/',addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
  df <- t(sapply(0:20, function(i){
    print(i)
    expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=i)
    r <- exprdiff(expr=expr2,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=FALSE)
    saveRDS(r, paste0('./addsignal/result/limma/',addSignalType,'/',geneProp,'_',i,'.rds'))
    r <- r[order(r[,3]),]
    sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
    res <- list()
    res[['res']] <- r
    res[['sensfdr']] <- c(i,AreaUnderSensFdr(sensfdr))
    print(c(i,AreaUnderSensFdr(sensfdr)))
    return(c(i,AreaUnderSensFdr(sensfdr)))
  }))
  colnames(df)[1] <- c('Parameter')
  saveRDS(df, paste0('./addsignal/result/limma/',addSignalType,'/',geneProp,'_Para_FdrDiff_Area_Group1256.rds'))
} 

if (method == 'permu'){
  dir.create(paste0('./addsignal/result/permu/',addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
  expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)
  res <- exprdiff(expr=expr2,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE, permutime=10000)
  r <- res[['res']]
  r <- r[order(r[,3]),]
  sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
  v <- c(addSignalPara,AreaUnderSensFdr(sensfdr))
  names(v)[1] = 'Parameter'
  res[['sensfdr']] <- v
  saveRDS(res, paste0('./addsignal/result/permu/', addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))
} 

if (method == 'permu_IR'){ ## ignore repetitive
  dir.create(paste0('./addsignal/result/permu_IR/',addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
  expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)
  res <- exprdiff(expr=expr2,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE, permutime=10000, IgnoreRepetitive = TRUE)
  r <- res[['res']]
  r <- r[order(r[,3]),]
  sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
  v <- c(addSignalPara,AreaUnderSensFdr(sensfdr))
  names(v)[1] = 'Parameter'
  res[['sensfdr']] <- v
  saveRDS(res, paste0('./addsignal/result/permu_IR/', addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))
} 

if (grepl('tradeSeq', method)){
  # dir.create(paste0('./addsignal/result/', method, '/',addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(paste0('./addsignal/result/', method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))){
    expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)
    counts <- round(exp(expr2 + 1))
  
    psn <- seq(1, length(pseudotime))
    names(psn) <- pseudotime
    pdt <- data.frame(curve1 = psn, curve2 = psn)
    rownames(pdt) <- names(psn)
    pdt = pdt[colnames(counts), ]
  
    v <- (sub(':.*','',allp) %in% paste0('BM',c(1,2,5,6)) + 0)
    v <- ifelse(v==1, 0.99, 0.01)
    cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
    rownames(cellWeights) <- colnames(counts)
  
    set.seed(12345)
    sce <- fitGAM(counts = counts, pseudotime = pdt, cellWeights = cellWeights,
                     nknots = 6, verbose = FALSE,parallel=TRUE)
    saveRDS(sce, paste0('./addsignal/result/', method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))
  } else {
    sce <- readRDS(paste0('./addsignal/result/', method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'_sce.rds'))
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
    saveRDS(Final, paste0('./addsignal/result/', method,'/',addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  
  }
    
}

rm(list=ls())
