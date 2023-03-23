library(rhdf5)
selid <- readRDS('/home/whou10/scratch16/whou10/trajectory_variability/tb/res/sex_batchcase/selid.rds')
pid <- 2
ddir = '/home/whou10/scratch16/whou10/trajectory_variability/tb/data/sex/'

pt <- readRDS(paste0(ddir,'ptpc',pid,'.rds'))
cellanno <- readRDS(paste0(ddir, 'cellanno.rds'))
design <- readRDS(paste0(ddir, 'design.rds'))

design <- design[sub('_.*','',rownames(design)) %in% selid,]
design <- design[,colSums(design) > 0]
tb <- read.table('/home/whou10/scratch16/whou10/trajectory_variability/tb/data/raw/GSE158769_meta_data.txt', sep = '\t',header=T)
tb <- unique(tb[,c('TB_status','donor')])
design <- cbind(tb=as.numeric(tb$TB_status[match(sub('_.*','',rownames(design)),tb$donor)]=='CASE'),design)

cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
pt <- pt[cellanno[,1]]

## hdf5
library(parallel)
source('/home/whou10/scratch16/whou10/trajectory_variability/h5func/01_function.R')
af <- list.files('/home/whou10/scratch16/whou10/trajectory_variability/h5func',pattern = 'multi')
for (f in af) source(paste0('/home/whou10/scratch16/whou10/trajectory_variability/h5func/',f))
res <- testpt(expr=paste0(ddir, 'exprpc',pid,'.h5'), cellanno=cellanno, pseudotime=pt, design=design,testvar=ncol(design), test.type = 'Variable', demean = FALSE, overall.only = F, test.method = 'permutation', ncores = 20)
saveRDS(res,file=paste0('/home/whou10/scratch16/whou10/trajectory_variability/tb/res/sex_batchcase/lamian_pm2.rds'))


