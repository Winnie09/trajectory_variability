library(here)
setwd(here())
source('function/01_function.R')
for (path in c('lymph', 'erythroid', 'monocyte')){
  ddir <- paste0('hca/real/build_from_tree_variability/result/', path, '/')
  # rdir <- paste0('hca/real/testvar/result/EM_pm/', path, '/')
  cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
  rownames(cellanno) = cellanno[,1]
  cellanno[,2] <- sub(':.*','',cellanno[,2])
  design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
  rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/result/', path, '/')
  dir.create(rdir, recursive = T)
  for (prop in c(0.01,0.05, 0.25,seq(0.1, 0.9, 0.1))) {
    print(paste0(path, '_', prop))
  pseudotime = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/data/',path, '/', prop, '.rds'))
  res <- cellPropTest(cellanno = cellanno[names(pseudotime), ], pseudotime = pseudotime, design = design, ncores = 10, test.type = 'variable')
  saveRDS(res, paste0(rdir, 'cell_proportion_test_', prop,'.rds'))  
  }
}

