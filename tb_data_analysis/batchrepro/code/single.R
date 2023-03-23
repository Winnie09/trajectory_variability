library(rhdf5)
library(parallel)
pid <- 2
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

seed <- as.numeric(commandArgs(trailingOnly = T))
print(seed)
set.seed(seed)
us <- unique(cellanno[,2])
batch1 <- grep('batch',colnames(design),value=T)[seed]
batch2 <- setdiff(grep('batch',colnames(design),value=T),batch1)[1:10]
samp1 <- grep(batch1,us,value=T)
samp2 <- unlist(sapply(batch2,grep,us,value=T))
cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
cellanno2 <- cellanno[cellanno[,2] %in% samp2,]

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/01_function.R')
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func',pattern = 'multi')
for (f in af) source(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/',f))

des <- cbind(1,design[samp1,ncol(design)])
colnames(des)[ncol(des)] <- 'sex'

res <- testpt(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/exprpc',pid,'.h5'), cellanno1, pt[cellanno1[,1]], design=des,testvar=ncol(des),test.type='Variable')

saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/batchrepro/res/single/',seed,'.rds'))




