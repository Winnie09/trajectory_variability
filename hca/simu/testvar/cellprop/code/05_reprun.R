library(here)
setwd(here())
source('function/01_function.R')
path <- 'erythroid'
ddir <- paste0('hca/real/build_from_tree_variability/result/', path, '/')
# rdir <- paste0('hca/real/testvar/result/EM_pm/', path, '/')
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
rownames(cellanno) = cellanno[,1]
cellanno[,2] <- sub(':.*','',cellanno[,2])
design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/represult/')
dir.create(rdir, recursive = T)
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/repdata/') 
prop <- commandArgs(trailingOnly = T)
af <- af[grep(paste0('_',prop,'.rds'),af)]
res <- mclapply(af,function(f) {
  print(f)
  pseudotime = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/repdata/',f))
  res <- cellPropTest(cellanno = cellanno[names(pseudotime), ], pseudotime = pseudotime, design = design, 
                      ncores = 10, test.type = 'variable')
  res[[1]]
},mc.cores=detectCores()) 
saveRDS(res, paste0(rdir, prop,'.rds'))  
  
