rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
dir.create(pdir)

## overall performance
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf.rds')
pd[,2] <- as.character(pd[,2])
pd <- pd[pd[,2]!='EM_chisq', ]
pd <- pd[pd[,2]!='TSCAN', ]
library(ggplot2)
library(gridExtra)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
scale_color_brewer(palette = 'Dark2') 
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('Area under sensitivity-realFDR curve') +
scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc.pdf'),width=7,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()


## mean only
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_meanOnly.rds')
pd[,2] <- as.character(pd[,2])
pd <- pd[pd[,2]!='EM_chisq', ]
pd <- pd[pd[,2]!='TSCAN', ]
str(pd)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('Area under sensitivity-realFDR curve') +
  scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc_meanOnly.pdf'),width=7,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

## trend Only
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendOnly.rds')
pd[,2] <- as.character(pd[,2])
pd <- pd[pd[,2]!='EM_chisq', ]
pd <- pd[pd[,2]!='TSCAN', ]
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  # ylab('Area Under Sensitivity-RealFDR Curve') +
    ylab('Area under sensitivity-realFDR curve') +
  scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc_trendOnly.pdf'),width=7,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

## trend Mean
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendMean.rds')
pd[,2] <- as.character(pd[,2])
pd <- pd[pd[,2]!='EM_chisq', ]
pd <- pd[pd[,2]!='TSCAN', ]
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('FDR.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
    ylab('Area under sensitivity-realFDR curve') +
  scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc_trendMean.pdf'),width=7,height=3)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()



pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf.rds')
pd[,2] <- as.character(pd[,2])
pd0 <- pd[pd[,2]!='EM_chisq', ]

pd0$comparison <- 'overall'
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_meanOnly.rds')
pd[,2] <- as.character(pd[,2])
pd1 <- pd[pd[,2]!='EM_chisq', ]

pd1$comparison = 'mean only'
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendOnly.rds')
pd[,2] <- as.character(pd[,2])
pd2 <- pd[pd[,2]!='EM_chisq', ]

pd2$comparison = 'trend only'
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendMean.rds')
pd[,2] <- as.character(pd[,2])
pd3 <- pd[pd[,2]!='EM_chisq', ]

pd3$comparison = 'trend & mean'
pd <- rbind(pd0, pd1,pd2,pd3)
pd <- pd[pd[,2]!='TSCAN', ]
pd$comparison <- factor(pd$comparison, levels = c('overall', 'trend only', 'trend & mean','mean only'))
pd[pd[,2] == 'EM_pm', 2] <- 'Lamian' ## method name

pdf(paste0(pdir, 'compare_fdr_diff_all.pdf'),width=5.5,height=2.8)
ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_bw() + 
  xlim(c(0,4))+
  ylab('FDR.difference') +
  scale_color_brewer(palette = 'Set1')  +
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
  scale_color_brewer(palette = 'Set1')  +
  facet_wrap(~comparison, nrow =1)+
    theme(legend.position = 'bottom')
dev.off()

