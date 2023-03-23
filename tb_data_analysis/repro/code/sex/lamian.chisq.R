library(rhdf5)
pid <- 2
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

seed <- as.numeric(commandArgs(trailingOnly = T))
set.seed(seed)
samp1 <- sample(rownames(design),nrow(design)/2)
samp2 <- setdiff(rownames(design),samp1)
cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
cellanno2 <- cellanno[cellanno[,2] %in% samp2,]

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/01_function.R')
res <- testpt(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/exprpc',pid,'.h5'), cellanno1, pt[cellanno1[,1]], design=design[samp1,],test.type='Variable',test.method = 'chisq')
saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/lamian_chisq/',seed,'_1.rds'))

res <- testpt(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/exprpc',pid,'.h5'), cellanno2, pt[cellanno2[,1]], design=design[samp2,],test.type='Variable',test.method = 'chisq')
saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/lamian_chisq/',seed,'_2.rds'))

