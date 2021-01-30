rm(list=ls())
geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
## be reminded that the pca is done on cv > 0.5 (for all samples) genes. You should redo it !!!
# geneProp <- 0.05
# addSignalType <- 'linear'
# addSignalPara <-  5

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
rdir <- './addsignal/data/plain/1256/'
if (!file.exists(paste0(rdir, addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))){
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
  sample <- sub(':.*','',colnames(expr))
  names(sample) <- colnames(expr)
  design <- cbind(c(1,1,0,0,1,1,0,0))
  row.names(design) <- paste0('BM',1:8)
  colnames(design) <- 'group'
  
  ## construct pseudotime 
  library(TSCAN)
  pseudotime <- TSCANorder(exprmclust(t(dr),reduce = F,clusternum=4,clustermethod='kmeans'),orderonly = T)
  psn <- 1:length(pseudotime)
  names(psn) <- pseudotime
  order = data.frame(Pseudotime = seq(1, length(pseudotime)), Cell=pseudotime)
  
  ## add signal
  set.seed(12345)
  selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))
  celln <- sapply(colnames(expr), function(i) paste0(strsplit(i,':')[[1]][2],':',strsplit(i,':')[[1]][3]))
  set.seed(12345)
  expr3 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)
    
  ## save
  dir.create(paste0(rdir, addSignalType), recursive = TRUE)
  saveRDS(list(expr = expr3, selgene = selgene, pseudotime = pseudotime, design = design), paste0(rdir, addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  
}

