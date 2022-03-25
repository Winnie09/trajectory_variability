library(Seurat)
pid = 2
pt = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
expr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds')
cellanno = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds') ## ## female = 1, male = 0

cc.s = colSums(expr[intersect(cc.genes[[1]], rownames(expr)),])
cc.g2m = colSums(expr[intersect(cc.genes[[2]], rownames(expr)),])

mat = rbind(cc.s, cc.g2m)
colnames(mat) = colnames(expr)

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res = testpt(expr = mat, design = design, cellanno = cellanno, pseudotime = pt, test.type = 'Variable', ncores = 1)

saveRDS(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/res/tb_sex_testres.rds')
