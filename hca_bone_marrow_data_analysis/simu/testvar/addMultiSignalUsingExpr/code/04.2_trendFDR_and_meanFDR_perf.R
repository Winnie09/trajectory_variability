rm(list=ls())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
library(here)
setwd(here())
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds') # trend diff only
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds') # mean diff only
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds') # trend and mean
selgene = selgene1
rmgene = c(selgene2, selgene3)
str(selgene)
str(rmgene)
source('function/01_function.R')
library(gridExtra)

perflist <- list()
m = 'EM_pm'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.trendDiff'], zscore = r$statistics[,'z.trendDiff'], stringsAsFactors = F)
  res = res[complete.cases(res),]
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'trend.test', AreaUnderSensFdr(SensFdr(c(selgene1,selgene3), res)))
})
perflist[['trend.test']] = t(df1.nc)

df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.meanDiff'], zscore = r$statistics[,'z.meanDiff'], stringsAsFactors = F)
  res = res[complete.cases(res),]
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'mean.test', AreaUnderSensFdr(SensFdr(c(selgene2,selgene3), res)))
})
# df1.nc[3,4] = (as.numeric(df1.nc[3,3]) + as.numeric(df1.nc[3,5]))/2
# df1.nc[4,4] = (as.numeric(df1.nc[4,3]) + as.numeric(df1.nc[4,5]))/2
perflist[['mean.test']] = t(df1.nc)

###
af = list.files(paste0(ddir, 'EM_chisq/'))
af <- af[grepl('.rds', af)]
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, 'EM_chisq/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.chisq.trendDiff'], stat = (r$ll3-r$ll2)/r$statistics[,'df.diff.trendDiff'], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'trend.test', AreaUnderSensFdr(SensFdr(c(selgene1, selgene3), res)), 'Lamian.chisq')
})
perflist[['Lamian.chisq.trend.test']] = t(df1)

df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, 'EM_chisq/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.chisq.meanDiff'], stat = (r$ll2-r$ll1)/r$statistics[,'df.diff.meanDiff'], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'mean.test', AreaUnderSensFdr(SensFdr(c(selgene2, selgene3), res)), 'Lamian.chisq')
})
perflist[['Lamian.chisq.mean.test']] = t(df1)



perf <- do.call(rbind,perflist)
perf <- data.frame('SignalStrength'=as.numeric(perf[,1]),'Method'=perf[,2],'Fdr.Diff'=as.numeric(perf[,3]),'AUC'=as.numeric(perf[,4]))
saveRDS(perf, paste0(rdir,'perf_trenddiff_meandiff.rds'))



