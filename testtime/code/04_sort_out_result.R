source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/selgene/selgene.rds')
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/result/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/result/perf/'
dir.create(rdir, showWarnings = F, recursive = T)

## EM_SelectKnots
m = 'EM_SelectKnots'
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('testres', af)]
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  c(sub('.rds','',f), 'EM_SelectKnots', AreaUnderSensFdr(SensFdr(selgene, r[['res']])))
})
df1 = t(df1)
  
## tradeSeq
m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('sce', af)]

df <- lapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  a = c(sub('.rds','',f), 'tradeSeq_startVsEndTest', AreaUnderSensFdr(SensFdr(selgene, r[['startVsEndTest']])))
  b = c(sub('.rds','',f), 'tradeSeq_associationTest', AreaUnderSensFdr(SensFdr(selgene, r[['associationTest']])))
  rbind(a, b)
})
df2 <- do.call(rbind, df)

## tscan
m = 'tscan'
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('testres', af)]
df3 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  c(sub('.rds','',f),  'tscan', AreaUnderSensFdr(SensFdr(selgene, r[['res']])))
})
df3 = t(df3)

## monocle2
m = 'monocle2'
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('testres', af)]
df4 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  c(sub('.rds','',f),  'monocle2', AreaUnderSensFdr(SensFdr(selgene, r[['res']])))
})
df4 = t(df4)

## monocle3
m = 'monocle3'
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('testres', af)]
df5 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  c(sub('.rds','',f),  'monocle3', AreaUnderSensFdr(SensFdr(selgene, r)))
})
df5 = t(df5)

## concatenate
res <- rbind(df1, df2, df3, df4, df5)
colnames(res) <- c('Type', 'Method', 'Fdr.Diff', 'AUC')
saveRDS(res, paste0(rdir, 'perf.rds'))
rm(list=ls())



