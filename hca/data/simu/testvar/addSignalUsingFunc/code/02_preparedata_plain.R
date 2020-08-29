geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
# geneProp <- 0.05
# addSignalType <- 'linear'
# addSignalPara <-  1

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testvar/data/plain/'
if (file.exists(paste0(rdir, addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))) break

source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')

expr <- readRDS('./testtime/data/data/null/hsc_mep_ery.rds')
sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'
pseudotime <- readRDS('./testtime/data/data/null/pseudotime.rds')

## add signal
set.seed(12345)
selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))
expr3 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara, type = 'Variable')

## save
dir.create(paste0(rdir, addSignalType), recursive = TRUE)
saveRDS(list(expr = expr3, selgene = selgene, design = design), paste0(rdir, addSignalType,'/', geneProp,'_', addSignalPara,'.rds'))  

