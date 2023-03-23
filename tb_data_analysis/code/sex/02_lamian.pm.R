library(rhdf5)
pid <- as.numeric(commandArgs(trailingOnly = T))
ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex_addBatch/'

pt <- readRDS(paste0(ddir,'ptpc',pid,'.rds'))
cellanno <- readRDS(paste0(ddir, 'cellanno.rds'))
design <- readRDS(paste0(ddir, 'design.rds'))

## hdf5
library(parallel)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/01_function.R')
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func',pattern = 'multi')
for (f in af) source(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/',f))
res <- testpt(expr=paste0(ddir, 'exprpc',pid,'.h5'), cellanno=cellanno, pseudotime=pt, design=design,testvar=ncol(design), test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'permutation', ncores = 5)
saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex_addBatch/pc',pid,'/lamian_pm.rds'))


