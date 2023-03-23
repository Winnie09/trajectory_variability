source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/saveh5.R')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/expr.rds')
expr <- expr[sample(rownames(expr),100),]
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
path <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/subexpr.h5'
saveh5(expr,pseudotime,cellanno,path)
