library(here)
setwd(here())
d  = readRDS('/home/whou10/scratch16/whou10/trajectory_variability/cellcycle/res/tb_cellcycle_score.rds')
cellanno = readRDS('tb/data/sex/cellanno.rds')

#library(Seurat)
pid = 2
pt = readRDS(paste0('tb/data/sex/ptpc',pid,'.rds'))
#expr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds')
#cellanno = readRDS('tb/data/sex/cellanno.rds')
design <- readRDS('tb/data/sex/design.rds') ## ## female = 1, male = 0

# cc.s = colSums(expr[intersect(cc.genes[[1]], rownames(expr)),])
# cc.g2m = colSums(expr[intersect(cc.genes[[2]], rownames(expr)),])

# mat = rbind(cc.s, cc.g2m)
# colnames(mat) = colnames(expr)

source('function/01_function.R')
res = testpt(expr = t(d[,1:3]), design = design, cellanno = cellanno, pseudotime = pt, test.type = 'Variable', ncores = 1)

saveRDS(res, 'cellcycle/res/tb_sex_testres.rds')


