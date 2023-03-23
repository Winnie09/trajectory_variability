

rm(list=ls())
library(here)
setwd(here())
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
m = 'EM_pm'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.overall'], zscore = r$statistics[,'z.overall'], stringsAsFactors = F)
  res = res[order(res[,1], res[,2]), ]
  
  res = res[!rownames(res) %in% rmgene, ]
  c(sub('.rds','',f), 'EM_pm', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_pm']] = t(df1.nc)


m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))[[1]]
  res = res[order(res[, 3],-abs(res[, 1])),]
  c(sub('.rds', '', f), 'tradeSeq_diffEndTest', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tradeSeq_diffEndTest']] = t(df3)

m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))[[2]]
  res[is.na(res[,3]),3] <- 1
  res[is.na(res[,2]),2] <- 1
  res[is.na(res[,1]),1] <- 0
  res = res[order(res[, 3],-abs(res[, 1])),]
  c(sub('.rds', '', f), 'tradeSeq_patternTest', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tradeSeq_patternTest']] = t(df3)

m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))[[3]]
  res[is.na(res[,3]),3] <- 1
  res[is.na(res[,2]),2] <- 1
  res[is.na(res[,1]),1] <- 0
  res = res[order(res[, 3],-abs(res[, 1])),]
  c(sub('.rds', '', f), 'tradeSeq_earlyDETest', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tradeSeq_earlyDETest']] = t(df3)



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


saveRDS(perflist, paste0(rdir,'perflist_trendMean.rds'))
perf <- do.call(rbind, perflist)
colnames(perf) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)

perf[,1] <- as.numeric(as.character(perf[,1]))
perf[,3] <- as.numeric(as.character(perf[,3]))
perf[,4] <- as.numeric(as.character(perf[,4]))

saveRDS(perf, paste0(rdir,'perf_trendMean.rds'))
rm(list=ls())








