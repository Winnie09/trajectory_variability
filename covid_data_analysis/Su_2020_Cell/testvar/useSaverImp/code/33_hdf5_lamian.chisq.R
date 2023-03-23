library(rhdf5)
library(here)
setwd(here())
cellanno = readRDS('covid/Su_2020_Cell/data/cellanno.rds')
pseudotime = readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
design = readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds')

path <- 'covid/Su_2020_Cell/testvar/useSaverImp/data/expr.h5'
source('function/01_function.R')
source('h5func/01_function.R')
res <- testpt(path, cellanno, pseudotime, design=design,test.type='Variable',test.method='chisq')
saveRDS(res,file='covid/Su_2020_Cell/testvar/useSaverImp/result/lamian.chisq.hdf5/Mod_Mi/res.rds')


