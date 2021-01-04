library(here)
setwd(here())
ddir <- 'hca/simu/testvar/nullsimu/result/'
rdir <- 'hca/simu/testvar/nullsimu/result/perf/'
dir.create(rdir)

res.tscan <- readRDS(paste0(rdir, 'tscan/res.rds'))
res.em <- readRDS(paste0(rdir, 'EM_NOT_centered/res.rds'))
res.mean <- readRDS(paste0(rdir, 'meandiff/res.rds'))


num <- c(sum(res.tscan[,3]<0.05), sum(res.mean[,5]<0.05), 0)
names(num) <- c('TSCAN', 'Limma', 'EM_NOT_centered')
saveRDS(num, paste0(rdir, 'number_of_false_positives.rds'))
