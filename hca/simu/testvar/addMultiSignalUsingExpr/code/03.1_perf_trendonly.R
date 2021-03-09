rm(list=ls())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
library(here)
setwd(here())
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds')
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds')
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds')
selgene = selgene1
rmgene = c(selgene2, selgene3)
str(selgene)
str(rmgene)
source('function/01_function.R')

perflist <- list()
m = 'EM_pm'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = r$statistics
  res = res[order(res[,1], -res[,2]), c(1:3)]
  colnames(res) <- c('fdr', 'zscore', 'pvalue')
  res = res[!rownames(res) %in% rmgene, ]
  c(sub('.rds','',f), 'EM_pm', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_pm']] = t(df1.nc)


m = 'EM_chisq'
af = list.files(paste0(ddir, m, '/'))
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = r$statistics
  res = res[order(res[,1], res[,2]), c(1:2)]
  colnames(res) <- c('fdr', 'pvalue')
  res = res[!rownames(res) %in% rmgene, ]
  c(sub('.rds','',f), 'EM_chisq', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_chisq']] = t(df1)


m = 'meandiff'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[!rownames(res) %in% rmgene, ]
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'limma', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['limma']] = t(df2)


m = 'tscan'
af = list.files(paste0(ddir, m, '/'))
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[!rownames(res) %in% rmgene, ]
  res = res[order(res[, 3],-abs(res[, 2])),]
  c(sub('.rds', '', f), 'TSCAN', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tscan']] = t(df3)

## concatenate
saveRDS(perflist, paste0(rdir,'perflist_trendOnly.rds'))
perf <- do.call(rbind, perflist)
colnames(perf) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)

perf[,1] <- as.numeric(as.character(perf[,1]))
perf[,3] <- as.numeric(as.character(perf[,3]))
perf[,4] <- as.numeric(as.character(perf[,4]))

saveRDS(perf, paste0(rdir,'perf_trendOnly.rds'))
rm(list=ls())


