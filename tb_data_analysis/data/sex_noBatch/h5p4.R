source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/saveh5.R')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc4.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
path <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/exprpc4.h5'
saveh5(expr,pseudotime,cellanno,path)

