library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/saveh5.R')
ddir <- 'hca/real/build_from_tree_variability/result/monocyte/'
rdir <- 'hca/real/testvar/data/data/monocyte/'
dir.create(rdir, recursive = T)
expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
path <- paste0(rdir, 'expr.h5')
saveh5(expr,pseudotime,cellanno,path)


library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/saveh5.R')
ddir <- 'hca/real/build_from_tree_variability/result/erythroid/'
rdir <- 'hca/real/testvar/data/data/erythroid/'
dir.create(rdir, recursive = T)
expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
path <- paste0(rdir, 'expr.h5')
saveh5(expr,pseudotime,cellanno,path)


library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/h5func/saveh5.R')
ddir <- 'hca/real/build_from_tree_variability/result/lymph/'
rdir <- 'hca/real/testvar/data/data/lymph/'
dir.create(rdir, recursive = T)
expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
path <- paste0(rdir, 'expr.h5')
saveh5(expr,pseudotime,cellanno,path)



