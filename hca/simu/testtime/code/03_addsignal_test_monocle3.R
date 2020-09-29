# --------------
# global setting
# --------------
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
ddir <- './hca/data/simu/testtime/'
rdir <- './hca/simu/testtime/result/addsignal/'
method = 'monocle3'
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))

# ------------
# prepare data
# ------------
pt <- readRDS('./hca/data/simu/testtime/null/pseudotime.rds')
pseudotime <- pt[,2]
names(pseudotime) <- pt[,1]

# ----------------------------------
# test  [one group along pseudotime]
# ----------------------------------
for (clusterType in 1:4){
  print(clusterType)
  for (pctGene in 1:4){
    print(pctGene)
    library(spdep)
    expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
    expr <- expr[, names(pseudotime)]
    pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/null/hsc_mep_ery_integrated_pca.rds')
    res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
    saveRDS(res, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
  }
}
  
    
for (clusterType in 5:9){
  print(clusterType)
  for (pctGene in 1:4){
    print(pctGene)
    library(spdep)
    expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
    expr <- expr[, names(pseudotime)]
    pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/null/hsc_mep_ery_integrated_pca.rds')
    res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
    saveRDS(res, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
  }
}
  

