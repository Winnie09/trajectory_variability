library(rhdf5)
library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/lymph/'
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)

path <- 'hca/real/testvar/data/data/monocyte/expr.h5'
source('function/01_function.R')
source('h5func/01_function.R')
res <- testpt(path, cellanno, pseudotime, design=design,test.type='Variable')
saveRDS(res,file='hca/real/testvar/result/lamian.pm.hdf5/monocyte/gender/res.rds')

