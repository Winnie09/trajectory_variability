######
rm(list=ls())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
library(here)
setwd(here())
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/perf/'
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/data/selgene/selgene1.rds')
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/data/selgene/selgene2.rds')
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/data/selgene/selgene3.rds')
selgene = selgene3
rmgene = c(selgene1, selgene2)

str(selgene)
str(rmgene)
source('function/01_function.R')

perflist <- list()
m = 'Lamian.pm'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.overall'], zscore = r$statistics[,'z.overall'], stringsAsFactors = F)
  res = res[order(res[,1], res[,2]), ]
  
  res = res[!rownames(res) %in% rmgene, ]
  c(sub('.rds','',f), 'Lamian', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['Lamian']] = t(df1.nc)

m = 'limma'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[!rownames(res) %in% rmgene, ]
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'limma', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['limma']] = t(df2)

m = 'Lamian.chisq'
af = list.files(paste0(ddir, m, '/'))
af <- af[grepl('.rds', af)]
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  if ('df.diff.overall' %in% colnames(r$statistics)){
    res = data.frame(fdr = r$statistics[,'fdr.chisq.overall'], stat = (r$ll3-r$ll1)/r$statistics[,'df.diff.overall'], stringsAsFactors = F)
    res = res[order(res[,1], -res[,2]), ]
  } else {
    res = data.frame(fdr = r$statistics[,'fdr.chisq.overall'], stat = (r$ll3-r$ll1)/r$statistics[,'pval.chisq.overall'], stringsAsFactors = F)
    res = res[order(res[,1], res[,2]), ]
  }

  res = res[!rownames(res) %in% rmgene, ]
  c(sub('.rds','',f), 'Lamian.chisq', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['Lamian.chisq']] = t(df1)

## concatenate
saveRDS(perflist, paste0(rdir,'perflist_trendMean.rds'))
perf <- do.call(rbind, perflist)
colnames(perf) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)
perf[,1] <- as.numeric(as.character(perf[,1]))
perf[,3] <- as.numeric(as.character(perf[,3]))
perf[,4] <- as.numeric(as.character(perf[,4]))
saveRDS(perf, paste0(rdir,'perf_trendMean.rds'))



