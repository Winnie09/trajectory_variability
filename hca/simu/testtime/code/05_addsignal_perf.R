source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/selgene/selgene.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/addsignal/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/plot/addsignal/'
# perf = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/old_without_demean/result/addsignal/perf/perf.rds')
# perf[perf[,2]=='EM_SelectKnots',2] <- 'EM_NOT_centered'

## EM_centered, centered
af = list.files(paste0(rdir, 'EM_centered/'))
f = af[1]
df1<- t(sapply(af, function(f){
  print(f)
  Res = readRDS(paste0(rdir, 'EM_centered/', f))
  res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
}))

## EM_NOT_centered, centered
rdir2 <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/old_without_demean/result/addsignal/EM_SelectKnots/'
af = list.files(rdir2)
df6<- t(sapply(af, function(f){
  print(f)
  Res = readRDS(paste0(rdir2, f))
  res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_NOT_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
}))

## tradeSeq
m = 'tradeSeq'
af = list.files(paste0(rdir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('sce', af)]

df <- lapply(af, function(f){
  print(f)
  r = readRDS(paste0(rdir, m, '/', f))
  a = c(sub('.rds','',f), 'tradeSeq_startVsEndTest', AreaUnderSensFdr(SensFdr(selgene, r[['startVsEndTest']])))
  b = c(sub('.rds','',f), 'tradeSeq_associationTest', AreaUnderSensFdr(SensFdr(selgene, r[['associationTest']])))
  rbind(a, b)
})
df2 <- do.call(rbind, df)

## tscan
m = 'tscan'
af = list.files(paste0(rdir, m, '/'))
df3 <- t(sapply(af, function(f){
  print(f)
  r = readRDS(paste0(rdir, m, '/', f))
  c(sub('.rds','',f),  'tscan', AreaUnderSensFdr(SensFdr(selgene, r)))
}))


## monocle2
m = 'monocle2'
af = list.files(paste0(rdir, m, '/'))
df4 <- t(sapply(af, function(f){
  print(f)
  r = readRDS(paste0(rdir, m, '/', f))
  c(sub('.rds','',f),  'monocle2', AreaUnderSensFdr(SensFdr(selgene, r)))
}))


## monocle3
m = 'monocle3'
af = list.files(paste0(rdir, m, '/'))
df5 <- t(sapply(af, function(f){
  print(f)
  r = readRDS(paste0(rdir, m, '/', f))
  c(sub('.rds','',f),  'monocle3', AreaUnderSensFdr(SensFdr(selgene, r)))
}))


## concatenate
perf <- rbind(df1, df2, df3, df4, df5, df6)
colnames(perf) <- c('Type', 'Method', 'Fdr.Diff', 'AUC')
saveRDS(perf, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/addsignal/perf/perf.rds')
rm(list=ls())

# ----------------
# plot performance
# ----------------
pd <- data.frame(perf, SignalType = gsub('_.*','',sub('clusterType','',perf[,1])), stringsAsFactors = F)
pd$SignalStreghth <- as.numeric(sapply(pd$Type, function(i) sub('.*_', '', i) ))
pd[,2] = as.factor(pd[,2])
pd[,5] = as.factor(pd[,5])
pd[,3] = as.numeric(pd[,3])
pd[,4] = as.numeric(pd[,4])
library(ggplot2)
library(gridExtra)
pdf(paste0(pdir, 'compare_fdr_diff.pdf'),width=8,height=4)
ggplot(pd, aes(x = SignalStreghth, y = Fdr.Diff, color=Method)) + 
  geom_point(size=1)  + 
  geom_line(size=0.1) + 
  theme_classic() + 
  ylab('fdr.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Dark2') + 
  facet_wrap(~SignalType, scales = 'free')
dev.off()
pdf(paste0(pdir, 'compare_auc.pdf'), width=8,height=4)
ggplot(pd, aes(x = SignalStreghth, y = AUC, color=Method)) + 
  geom_point(size=1)  + 
  geom_line(size=0.1) + 
  theme_classic() + 
  ylab('Area Under Sensitivity Real_FDR curve') +
  scale_color_brewer(palette = 'Dark2') + 
  facet_wrap(~SignalType, scales = 'free')
dev.off()



