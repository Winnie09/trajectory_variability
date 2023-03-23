library(here)
setwd(here())
library(rhdf5)
pid <- 2
ddir = 'tb/data/sex/'
pt <- readRDS(paste0(ddir,'ptpc',pid,'.rds'))
cellanno <- readRDS(paste0(ddir, 'cellanno.rds'))
design <- readRDS(paste0(ddir, 'design.rds'))
# design = design[, c(1, 39, 2:38)]
design = design[,1, drop = F]

source('function/01_function.R')
af = list.files('function', pattern = 'multi.R$')
for (f in af){
  source(paste0('function/', f))
}
res <- cellPropTest(cellanno=cellanno, pseudotime=pt, design=design, test.type = 'time')

dir.create('tb/tcd/res/', recursive = T, showWarnings = F)
saveRDS(res,file='tb/tcd/res/pc2_res.rds')

