library(here)
setwd(here())
library(parallel)
source('function/01_function.R')
path <- 'erythroid'
# rdir <- paste0('hca/real/testvar/result/EM_pm/', path, '/')
cellanno = readRDS(paste0('hca/real/build_from_tree_variability/result/', path, '/input_cellanno.rds'))
rownames(cellanno) = cellanno[,1]
cellanno[,2] <- sub(':.*','',cellanno[,2])
design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/cellprop/result/testres/')
dir.create(rdir, recursive = T)
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/cellprop/data/'
af <- list.files(ddir) 
prop <- commandArgs(trailingOnly = T)[[1]]  ########
print(paste0('prop = ', prop))
af <- af[grep(paste0('_',prop,'.rds'),af)]
r <- lapply(af,function(f) {
  print(f)
  pseudotime = readRDS(paste0(ddir,f))
  res <- cellPropTest(cellanno = cellanno[names(pseudotime), ], pseudotime = pseudotime, design = design, 
                      ncores = 4, test.type = 'time')
  res[[1]]
}) 
print(str(r))
saveRDS(r, paste0(rdir, prop,'.rds'))  
  


