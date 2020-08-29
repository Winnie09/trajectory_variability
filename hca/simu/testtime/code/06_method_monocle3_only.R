setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('./function/01_function.R')
datadir <- './testtime/data/data/'
rdir <- './testtime/result/'
ddir <- './testtime/data/data/'
method <- 'monocle3'
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
dir.create(paste0(rdir, method), showWarnings = FALSE, recursive = TRUE)
pseudotime <- readRDS('./testtime/data/data/null/pseudotime.rds')
## one group along pseudotime

for (clusterType in seq(1, 10)){
  for (pctGene in seq(1, 4)){
    library(spdep)
    expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))
    psn <- as.numeric(pseudotime[,2])
    names(psn) <- pseudotime[,1]
    expr <- expr[, pseudotime[,2]]
    pca = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/hsc_mep_ery_integrated_pca.rds')
    res <- monocle3_time(expr=expr, cell_coords = pca[,1:4])
    saveRDS(res, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  
  }
}

  


