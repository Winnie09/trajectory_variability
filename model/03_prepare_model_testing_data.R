##### E-M algorithm
rm(list=ls())
geneProp <- 0.05
addSignalType <- 'linear'
addSignalPara <-  2

# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')

## source functions
# source('./function/01_function.R')
allf <- list.files('./function')
allf = allf[!grepl('01_function.R',allf)]
for (f in allf){  
  source(paste0('./function/', f))
}

# source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
library(parallel)
library(splines)
library(limma)
library(RColorBrewer)
suppressMessages(library(SingleCellExperiment))

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
set.seed(12345)
selgene <- sample(rownames(expr), round(geneProp*nrow(expr)))

sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'
  
expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)
saveRDS(list(expr=expr2, pseudotime=pseudotime,selgene=selgene, design=design),'/Users/wenpinhou/Dropbox/trajectory_variability/model/model_testing_data.rds')


celln <- sapply(colnames(expr), function(i) paste0(strsplit(i,':')[[1]][2],':',strsplit(i,':')[[1]][3]))
set.seed(12345)
newsample <- sample(sample)
namedf <- data.frame(oldname = colnames(expr), newname = paste0(newsample, ':',celln), stringsAsFactors = FALSE)
colnames(expr) <- namedf[,2]
pseudotime <- namedf[match(pseudotime, namedf$oldname), 'newname']
expr3 <- AddSignal(expr = expr, sample = newsample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)
saveRDS(list(expr=expr3, pseudotime=pseudotime,selgene=selgene, design=design),'/Users/wenpinhou/Dropbox/trajectory_variability/model/model_shuffle_sample_testing_data.rds')

