rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
dir.create(pdir)

## overall performance
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf.rds')
pd <- pd[complete.cases(pd), ]
pd[,2] <- as.character(pd[,2])
# pd <- pd[pd[,2]!='EM_chisq', ]
# pd <- pd[pd[,2]!='TSCAN', ]
# pd[grep('phenopath',pd[,2]), 2] <- 'phenopath'
# pd[grep('monocle2_trajtest',pd[,2]), 2] <- 'monocle2_trajTest'
# pd[grep('EM_pm',pd[,2]), 2] <- 'Lamian'
library(RColorBrewer)
colv = c(brewer.pal(9,'Set1')[c(1:2)], brewer.pal(8, 'Dark2')[1:8])
names(colv) <- c('Lamian', 'Lamian.chisq', setdiff(unique(pd[,2]), c('Lamian', 'Lamian.chisq')))

library(ggplot2)
library(gridExtra)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
scale_color_manual(values = colv)  + theme_compact() 
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  ylab('Area under sensitivity-realFDR curve') + 
scale_color_manual(values = colv)  + theme_compact() 
pdf(paste0(pdir, 'compare_fdr_diff_auc.pdf'),width=8,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()


## mean only
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_meanOnly.rds')
pd[,2] <- as.character(pd[,2])
str(pd)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
  scale_color_manual(values = colv)  + theme_compact() 
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('Area under sensitivity-realFDR curve') +
  scale_color_manual(values = colv)  + theme_compact() 
pdf(paste0(pdir, 'compare_fdr_diff_auc_meanOnly.pdf'),width=8,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

## trend Only
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendOnly.rds')
pd[,2] <- as.character(pd[,2])
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
  scale_color_manual(values = colv)  + theme_compact() 
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  # ylab('Area Under Sensitivity-RealFDR Curve') +
    ylab('Area under sensitivity-realFDR curve') +
  scale_color_manual(values = colv)  + theme_compact() 
pdf(paste0(pdir, 'compare_fdr_diff_auc_trendOnly.pdf'),width=8,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

## trend Mean
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendMean.rds')
pd[,2] <- as.character(pd[,2])
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
  scale_color_manual(values = colv)  + theme_compact() 
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
    ylab('Area under sensitivity-realFDR curve') +
  scale_color_manual(values = colv)  + theme_compact() 
pdf(paste0(pdir, 'compare_fdr_diff_auc_trendMean.pdf'),width=7,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()



pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf.rds')
pd <- pd[complete.cases(pd), ]
pd[,2] <- as.character(pd[,2])
pd0 <- pd
pd0$comparison <- 'overall'


pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_meanOnly.rds')
pd[,2] <- as.character(pd[,2])
pd1 <- pd
pd1$comparison = 'mean only'

pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendOnly.rds')
pd <- pd[complete.cases(pd),]
pd[,2] <- as.character(pd[,2])
pd2 <- pd
pd2$comparison = 'trend only'


pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendMean.rds')
pd <- pd[complete.cases(pd),]
pd[,2] <- as.character(pd[,2])
pd3 <- pd
pd3$comparison = 'trend & mean'

pd <- rbind(pd0, pd1,pd2,pd3)
pd <- pd[pd[,2]!='TSCAN', ]
pd$comparison <- factor(pd$comparison, levels = c('overall', 'trend only', 'trend & mean','mean only'))


pdf(paste0(pdir, 'compare_fdr_diff_all.pdf'),width=5.5,height=2.8)
ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_bw() + 
  xlim(c(0,4))+
  ylab('FDR.difference') +
  scale_color_manual(values = colv)  + theme_compact() +
  facet_wrap(~comparison, nrow =1) +
  theme(legend.position = 'bottom')
dev.off()

pdf(paste0(pdir, 'compare_auc_all.pdf'),width=5.5,height=2.8)
ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_bw() + 
    xlim(c(0,4))+
    ylab('AUC (sensitivity-realFDR)') +
  scale_color_manual(values = colv)  + theme_compact() +
  facet_wrap(~comparison, nrow =1)+
    theme(legend.position = 'bottom')
dev.off()





