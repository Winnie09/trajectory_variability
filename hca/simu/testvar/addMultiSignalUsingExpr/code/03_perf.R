ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/perf/'

selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testvar/addMultiSignalUsingExpr/selgene/selgene.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

m = 'EM'
af = list.files(paste0(ddir, m, '/'))
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_trenddiff', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df1 = t(df1)


m = 'EM_meandiff'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'EM_meandiff', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df2 = t(df2)

# ------------------
# other methods 
# ------------------
## tradeSeq
m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('sce', af)]

df <- lapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = r[['earlyDETest']]
  a = c(sub('.rds','',f), 'tradeSeq_earlyDETest', AreaUnderSensFdr(SensFdr(selgene, res[['res']])))
  res = r[['patternTest']]
  b = c(sub('.rds','',f), 'tradeSeq_patternTest', AreaUnderSensFdr(SensFdr(selgene, res[['res']])))
  res = r[['diffEndTest']]
  c = c(sub('.rds','',f), 'tradeSeq_diffEndTest', AreaUnderSensFdr(SensFdr(selgene, res[['res']])))
  rbind(a, b, c)
})
df3 <- do.call(rbind, df)

## monocle2
m = 'monocle2'
af = list.files(paste0(ddir, m, '/'))
df4 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  c(sub('.rds','',f),  'monocle2', AreaUnderSensFdr(SensFdr(selgene, r)))
})
df4 = t(df4)

## monocle3
m = 'monocle3'
af = list.files(paste0(ddir, m, '/'))
df5 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  r = r[order(r[,3], -abs(r[,1])), ]
  c(sub('.rds','',f),  'monocle3', AreaUnderSensFdr(SensFdr(selgene, r)))
})
df5 = t(df5)


## tscan
m = 'tscan'
af = list.files(paste0(ddir, m, '/'))
df6 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  c(sub('.rds','',f),  'tscan', AreaUnderSensFdr(SensFdr(selgene, r)))
})
df6 = t(df6)

res <- rbind(df1, df2, df3, df4, df5, df6)
colnames(res) <- c('Signal.Strength', 'Method', 'Fdr.Diff', 'AUC')
saveRDS(res, paste0(rdir, 'perf.rds'))


