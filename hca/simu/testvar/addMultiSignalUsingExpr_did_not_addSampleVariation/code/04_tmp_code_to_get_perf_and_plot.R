library(here)
setwd(here())
trendgene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds')
meangene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds')
meantrendgene <-readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds')
selgene <- c(trendgene, meangene, meantrendgene)



ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/'
source('function/01_function.R')
m = 'chisq'
af = list.files(paste0(ddir, m, '/'))
df1 <- sapply(af, function(f) {
  Res <- readRDS(paste0(ddir, m, '/', f))
  res = Res$statistics
  res = res[complete.cases(res), ]
  res = res[order(res[,1], res[,2]), ]
  c(sub('.rds', '', f), 'chisq', AreaUnderSensFdr(SensFdr(selgene, res)))
})

  
m = 'tscan'
af = list.files(paste0(ddir, m, '/'))
df2 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[, 3],-abs(res[, 2])),]
  c(sub('.rds', '', f), 'TSCAN', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist = list()
perflist[['chisq']] <- t(df1)
perflist[['tscan']] = t(df2)



## concatenate
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
saveRDS(perflist, paste0(rdir,'perflist.rds'))
perf <- do.call(rbind, perflist)
colnames(perf) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)
saveRDS(perf, paste0(rdir,'perf.rds'))


perf <- as.data.frame(perf)
pd <- perf


######## pd
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
str(pd)
pd$Method <- as.character(pd$Method)
pd <- pd[pd$Method != 'trenddiff', ]

pd$Method[as.character(pd$Method) == 'EM_NOT_centered'] <- 'TrajDiffTest'
pd$Method[pd$Method == 'meandiff'] <- 'MeanDiffTest'
pd[,1] <- as.numeric(pd[,1])
pd[,3] = as.numeric(as.character(pd[,3]))
pd[,4] = as.numeric(as.character(pd[,4]))
library(ggplot2)
library(gridExtra)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  # ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
  ylab('RealFDR - ReportedFDR') +
  scale_color_brewer(palette = 'Dark2')  +
  ylim(c(-0.5, 0.5))
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  # ylab('Area Under Sensitivity-RealFDR Curve') +
  ylab('AUC') + 
  scale_color_brewer(palette = 'Dark2') +
  ylim(c(0,1))
pdf(paste0(pdir, 'compare_fdr_diff_auc.pdf'),width=7,height=2.2)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()




