rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds')
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds')
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds')
selgene = selgene3
rmgene = c(selgene1, selgene2)
str(selgene)
str(rmgene)
source('function/01_function.R')

perflist <- list()
# m = 'EM_centered'
# af = list.files(paste0(ddir, m, '/'))
# f = af[1]
# df1 <- sapply(af, function(f){
#   print(f)
#   r = readRDS(paste0(ddir, m, '/', f))
#   res = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
#   res = res[!rownames(res) %in% rmgene, ]
#   res = res[order(res[,1], -res[,2]), ]
#   
#   c(sub('.rds','',f), 'trenddiff', AreaUnderSensFdr(SensFdr(selgene, res)))
# })
# perflist[['EM_centered']] = t(df1)

m = 'EM_NOT_centered'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
  res = res[!rownames(res) %in% rmgene, ]
  res = res[order(res[,1], -res[,2]), ]
  
  c(sub('.rds','',f), 'EM_NOT_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_not_centered']] = t(df1.nc)

m = 'meandiff'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[!rownames(res) %in% rmgene, ]
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'meandiff', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['meandiff']] = t(df2)

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
saveRDS(perflist, paste0(rdir,'perflist_trendMean.rds'))
perf <- do.call(rbind, perflist)
colnames(perf) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)
saveRDS(perf, paste0(rdir,'perf_trendMean.rds'))
rm(list=ls())





