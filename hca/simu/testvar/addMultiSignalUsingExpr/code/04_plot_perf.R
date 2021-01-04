rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
dir.create(pdir)
## overall performance
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf.rds')
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
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  # ylab('Area Under Sensitivity-RealFDR Curve') +
  ylab('AUC') + 
  scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc.pdf'),width=7,height=2.2)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

## mean only
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_meanOnly.rds')
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
  ylab('RealFDR - ReportedFDR') +
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
   ylab('AUC') + 
  scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc_meanOnly.pdf'),width=7,height=2.2)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

## trend Only
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendOnly.rds')
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
  ylab('RealFDR - ReportedFDR') +
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  # ylab('Area Under Sensitivity-RealFDR Curve') +
   ylab('AUC') + 
  scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc_trendOnly.pdf'),width=7,height=2.2)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

## trend Mean
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendMean.rds')
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
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('AUC') +
  scale_color_brewer(palette = 'Dark2') 
pdf(paste0(pdir, 'compare_fdr_diff_auc_trendMean.pdf'),width=7,height=2.2)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()



