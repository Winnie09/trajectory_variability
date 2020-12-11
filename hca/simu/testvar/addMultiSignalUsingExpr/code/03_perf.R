library(here)
ddir <- paste0(here('hca','simu', 'testvar', 'addMultiSignalUsingExpr', 'result'), '/')
rdir <- paste0(here('hca','simu', 'testvar', 'addMultiSignalUsingExpr', 'result','perf'), '/')

selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testvar/addMultiSignalUsingExpr/selgene/selgene.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

perflist <- list()
m = 'EM_centered'
af = list.files(paste0(ddir, m, '/'))
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'trenddiff', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_centered']] = t(df1)

m = 'EM_NOT_centered'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_NOT_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['EM_not_centered']] = t(df1.nc)

m = 'meandiff'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'meandiff', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['meandiff']] = t(df2)


## concatenate
perf <- do.call(rbind, perflist)
colnames(perf) <- c('Type', 'Method', 'Fdr.Diff', 'AUC')
saveRDS(perf, paste0(rdir,'perf.rds'))
rm(list=ls())


