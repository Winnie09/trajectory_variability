rm(list=ls())
library(here)
setwd(here())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
selgene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene.rds')
source('function/01_function.R')

perflist <- list()
m = 'EM_NOT_centered'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,1], zscore = r$statistics[,2], stringsAsFactors = F)
  res = res[order(res[,1], res[,2]), ]
  c(sub('.rds','',f), 'EM_pm', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_pm']] = t(df1.nc)

m = 'chisq'
af = list.files(paste0(ddir, m, '/'))
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = r$statistics
  colnames(res)[1:2] <- c('fdr', 'pval')
  res = res[order(res[,1], res[,2]), ]
  c(sub('.rds','',f), 'EM_chisq', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_chisq']] = t(df1)

m = 'meandiff'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'limma', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['limma']] = t(df2)

m = 'tscan'
af = list.files(paste0(ddir, m, '/'))
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[, 3],-abs(res[, 2])),]
  c(sub('.rds', '', f), 'TSCAN', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tscan']] = t(df3)

## concatenate
saveRDS(perflist, paste0(rdir,'perflist.rds'))
perf <- do.call(rbind, perflist)
colnames(perf) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)
perf[,1] <- as.numeric(as.character(perf[,1]))
perf[,3] <- as.numeric(as.character(perf[,3]))
perf[,4] <- as.numeric(as.character(perf[,4]))
saveRDS(perf, paste0(rdir,'perf.rds'))
rm(list=ls())


