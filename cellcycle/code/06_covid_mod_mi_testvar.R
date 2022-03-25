library(Seurat)
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds') ## Mod = 0, Mi = 1

cc.s = colSums(expr[intersect(cc.genes[[1]], rownames(expr)),])
cc.g2m = colSums(expr[intersect(cc.genes[[2]], rownames(expr)),])

mat = rbind(cc.s, cc.g2m)
colnames(mat) = colnames(expr)

## subset mod and mi
cellanno = cellanno[cellanno[,2] %in% rownames(design), ]
selcell = cellanno[, 1]
mat = mat[, selcell]
pt = pt[selcell]

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res = testpt(expr = mat, design = design, cellanno = cellanno, pseudotime = pt, test.type = 'Variable', ncores = 1)

saveRDS(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/res/covid_Mod_Mi_testres.rds')

