# --------------
# global setting
# --------------
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
here()
method <- 'monocle3'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))
ddir <- here('hca/data/simu/testtime/addMultiSignalUsingExpr/')
rdir <- here('hca/simu/testtime/result/addsignal/')
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)

# ------------
# prepare data
# ------------
pt <- readRDS(here('hca/data/simu/testtime/poolSampleSignal/null/pseudotime.rds'))
pseudotime <- pt[,2]
names(pseudotime) <- pt[,1]

# ----------------------------------
# test  [one group along pseudotime]
# ----------------------------------
for (dataType in seq(1, 4)){
  library(spdep)
  expr <- readRDS(paste0(ddir, 'saver/', dataType, '.rds'))
  expr <- expr[, names(pseudotime)]
  pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/poolSampleSignal/null/hsc_mep_ery_integrated_pca.rds')
  res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
  saveRDS(res, paste0(rdir, method,'/', dataType,'.rds'))  
}

