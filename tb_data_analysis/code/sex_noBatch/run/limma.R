library(Matrix)
source('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/function/01_function.R')

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc2.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds')


res <- meandiff(expr = expr, cellanno = cellanno, design = design)
res <- res[order(res[,'adj.P.Val'], -abs(res[, 'logFC'])), ]
print(sum(res[,'adj.P.Val']<0.05))
saveRDS(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex/pc2/limma.rds')



