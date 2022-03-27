source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/saveh5.R')
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex_addBatch/'
expr <- readRDS(paste0(rdir, 'expr.rds'))
pseudotime <- readRDS(paste0(rdir, 'ptpc2.rds'))
cellanno <- readRDS(paste0(rdir, 'cellanno.rds'))
path <- paste0(rdir, 'exprpc2.h5')
saveh5(expr,pseudotime,cellanno,path)

