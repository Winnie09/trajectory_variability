library(Matrix)
source('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/function/01_function.R')

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/expr.rds')
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/case_design.rds')

res <- meandiff(expr = expr, cellanno = cellanno, design = design)
res <- res[order(res[,'adj.P.Val'], -abs(res[, 'logFC'])), ]
print(sum(res[,'adj.P.Val']<0.05))
saveRDS(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/case/limma.rds')


