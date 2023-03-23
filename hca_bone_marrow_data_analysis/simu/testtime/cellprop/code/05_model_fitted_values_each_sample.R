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

for (prop in seq(0, 0.5, 0.05)){
  af <- list.files(ddir) 
  print(paste0('prop = ', prop))
  af <- af[grep(paste0('_',prop,'.rds'),af)]
  f = af[1]
  print(f)
  pseudotime = readRDS(paste0(ddir,f))
  res <- cellPropTest(cellanno = cellanno[names(pseudotime), ], pseudotime = pseudotime, design = design, 
                      ncores = 4, test.type = 'time')
  pdir <- 'hca/simu/testtime/cellprop/plot/'
  dir.create(pdir, showWarnings = F, recursive = T)
  pdf(paste0(pdir, sub('.*_','',sub('.rds','',f)),'.pdf'), width = 10, height = 1.6)
  plotGene(res, 'prop',cellProp = T, plot.point = T, facet.sample = T, ncol = 8, point.size = 0.5, continuous = T)
  dev.off()  
}

