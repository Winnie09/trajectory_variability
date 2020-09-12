selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/selgene/selgene.rds')

perf = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/old_without_demean/result/addsignal/perf/perf.rds')
perf[perf[,2]=='EM_SelectKnots',2] <- 'EM_NOT_centered'

rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/addsignal/EM_SelectKnots/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/plot/addsignal/'

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

af = list.files(rdir)
f = af[1]
perf2 <- sapply(af, function(f){
  print(f)
  Res = readRDS(paste0(rdir, f))
  res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'EM_centered', AreaUnderSensFdr(SensFdr(selgene, res)))
})

perf <- rbind(perf, perf2)  
saveRDS(perf, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/addsignal/perf/perf.rds')

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


