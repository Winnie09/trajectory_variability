library(rhdf5)
pid <- as.numeric(commandArgs(trailingOnly = T))
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/01_function.R')
res <- testpt(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/exprpc',pid,'.h5'), cellanno, pt, design=design,test.type='Variable')
saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex/pc',pid,'/lamian_pm.rds'))

