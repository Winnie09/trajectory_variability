library(Matrix)
source('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/function/01_function.R')

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, cellanno[,1]]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

res <- meandiff(expr = expr, cellanno = cellanno, design = design)
res <- res[order(res[,'adj.P.Val'], -abs(res[, 'logFC'])), ]
print(sum(res[,'adj.P.Val']<0.05))
saveRDS(res, '/home-4/zji4@jhu.edu/scratch/diffpt/su/compres/limma/Mod_Mi.rds')

