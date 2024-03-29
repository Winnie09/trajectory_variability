library(here)
ddir <- paste0(here('hca','simu', 'testvar', 'addMultiSignalUsingExpr', 'result'), '/')
rdir <- paste0(here('hca','simu', 'testvar', 'addMultiSignalUsingExpr', 'perf'), '/')

selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testvar/addMultiSignalUsingExpr/selgene/selgene.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

m = 'EM_centered'
af = list.files(paste0(ddir, m, '/'))
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_trenddiff', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df1 = t(df1)

m = 'EM_NOT_centered'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_NOT_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df1.nc = t(df1.nc)

m = 'EM_meandiff'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'EM_meandiff', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df2 = t(df2)

## merge trend diff (centered) and mean diff, get the minimum fdr as the fdr for each gene
af = list.files(paste0(ddir, 'EM_centered/'))
df.merge <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, 'EM_centered/', f))
  res1 = data.frame(adj.P.Val = r[['fdr']], foldchange = r[['foldchange']], stringsAsFactors = F)
  r = readRDS(paste0(ddir, 'EM_meandiff/', f))
  res2 = r[,c('adj.P.Val', 'logFC')]
  colnames(res2) <- colnames(res1)
  res2 <- res2[rownames(res1), ]
  res <- lapply(rownames(res1), function(i){
    if (res1[i,1] <= res2[i,1]) {
      unlist(res1[i,,drop=T])
    } else {
      c(res2[i,1], res1[i,2])
    }
  })
  res <- do.call(rbind, res)
  rownames(res) <- rownames(res1)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_merge_min', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df.merge.min = t(df.merge)

## merge trend diff (centered) and mean diff, 1-(1-fdr1)(1-fdr2): under null,  the probability of oberving the data or more extreme cases in at least one situation
af = list.files(paste0(ddir, 'EM_centered/'))
df.merge <- sapply(af, function(f){
  print(f)
  r1 = readRDS(paste0(ddir, 'EM_centered/', f))
  
  r2 = readRDS(paste0(ddir, 'EM_meandiff/', f))
  r2 <- r2[names(r1[['fdr']]),]
  res <- data.frame(fdr = 1 - (1 - r1[['fdr']]) * (1 - r2[,'adj.P.Val']), 
                    foldchange = r1[['foldchange']],
                    stringsAsFactors = FALSE)
  rownames(res) <- names(r1[['fdr']])
  res = res[order(res[,1], res[,2]), ,drop = FALSE]
  c(sub('.rds','',f), 'EM_merge_atleast', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df.merge.atleast = t(df.merge)


## merge trend diff (centered) and mean diff, 1- fdr1 * fdr2: under null,  the probability of oberving the data or more extreme cases in at least one situation
af = list.files(paste0(ddir, 'EM_centered/'))
df.merge2 <- sapply(af, function(f){
  print(f)
  r1 = readRDS(paste0(ddir, 'EM_centered/', f))
  
  r2 = readRDS(paste0(ddir, 'EM_meandiff/', f))
  identical(names(r1[['fdr']]), rownames(r2[names(r1[['fdr']]),]))
  
  res <- data.frame(fdr = r1[['fdr']] * r2[names(r1[['fdr']]),'adj.P.Val'] -1, 
                    foldchange = r1[['foldchange']],
                    stringsAsFactors = FALSE)
  rownames(res) <- names(r1[['fdr']])
  res = res[order(res[,1], -res[,2]), ,drop = FALSE]
  c(sub('.rds','',f), 'EM_merge_atleast_1_fdr1fdr2', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df.merge.atleast_1_fdr1fdr2 = t(df.merge2)


## merge NOT-centered and mean diff, 1- fdr1 * fdr2: under null,  the probability of oberving the data or more extreme cases in at least one situation
af = list.files(paste0(ddir, 'EM_NOT_centered/'))
df.merge <- sapply(af, function(f){
  print(f)
  r1 = readRDS(paste0(ddir, 'EM_NOT_centered/', f))
  
  r2 = readRDS(paste0(ddir, 'EM_meandiff/', f))
  identical(names(r1[['fdr']]), rownames(r2[names(r1[['fdr']]),]))
  
  res <- data.frame(fdr = r1[['fdr']] * r2[names(r1[['fdr']]),'adj.P.Val'] -1, 
                    foldchange = r1[['foldchange']],
                    stringsAsFactors = FALSE)
  rownames(res) <- names(r1[['fdr']])
  res = res[order(res[,1], -res[,2]), ,drop = FALSE]
  # plot(res[,1], res[,2], xlab='1-fdr1*fdr2',ylab='trenddiff fc (orill-perll)', pch=20)
  c(sub('.rds','',f), 'EM_NOT_centered_1_fdr1fdr2', AreaUnderSensFdr(SensFdr(selgene, res)))
})
df.nc.1_fdr1fdr2 = t(df.merge)
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

## gather all results
res <- rbind(df1, df1.nc, df2, df3, df4, df5, df6)
colnames(res) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')

## add EM merge results
a = rbind(df.merge.min, df.merge.atleast, df.merge.atleast_1_fdr1fdr2, df.nc.1_fdr1fdr2)
b = data.frame(a[,1], a[,2], a[,3], a[,4], stringsAsFactors = F)
colnames(b) <- colnames(res)
Res = rbind(res,b)
saveRDS(Res, paste0(rdir, 'perf.rds'))
