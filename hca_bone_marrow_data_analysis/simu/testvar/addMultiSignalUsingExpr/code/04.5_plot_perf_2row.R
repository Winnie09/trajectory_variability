rm(list=ls())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
library(here)
setwd(here())
library(ggplot2)
library(RColorBrewer)
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
perf <- readRDS(paste0(rdir,'perf_trenddiff_meandiff.rds'))
pdonly <- perf

pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf.rds')
pd[,2] <- as.character(pd[,2])
pd[pd[,2]=='Lamian_chisq',2] <- 'Lamian.chisq'
pd0 <- pd
pd0$Comparison <- 'overall.test'

pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_meanOnly.rds')
pd[,2] <- as.character(pd[,2])
pd1 <- pd
pd1$Comparison = 'mean only'

pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendOnly.rds')
pd[,2] <- as.character(pd[,2])
pd[pd[,2]=='Lamian_chisq',2] <- 'Lamian.chisq'
pd2 <- pd
pd2$Comparison = 'trend only'


pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendMean.rds')
pd[,2] <- as.character(pd[,2])
pd3 <- pd
pd3$Comparison = 'trend & mean'

pd <- rbind(pd0, pd1,pd2,pd3)
pd$Comparison <- factor(pd$Comparison, levels = c('overall.test', 'trend only', 'trend & mean','mean only'))


pd <- rbind(pd,pdonly)
pd <- pd[complete.cases(pd), ]

pd[pd[,2]=='EM_pm',2] <- 'Lamian.pm'
pd[pd[,2]=='Lamian',2] <- 'Lamian.pm'
pd[,5] = factor(as.character(pd[,5]), levels = c("overall.test","trend.test", "mean.test", "trend only","trend & mean", "mean only"))


colv = c(brewer.pal(9,'Set1')[c(1:2)], brewer.pal(8, 'Dark2')[1:8])
names(colv) <- c('Lamian.pm', 'Lamian.chisq', setdiff(sort(unique(pd[,2])), c('Lamian.pm', 'Lamian.chisq')))

## save plot
pdf(paste0(pdir, 'compare_fdr_diff_auc_trenddiff_meandiff_2row_20220717.pdf'),width=7.5,height=4.3)
p1 <- ggplot(pd[pd[,2]!='monocle2_trajTest',], aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=0.5)  + 
  geom_line(size=0.3) + 
  xlim(c(0,4))+
  ylab('FDR.difference') +
  scale_color_manual(values = colv)  + 
  theme_classic() +
  theme(strip.text.x = element_text(size = 10), axis.title = element_text(size = 10), axis.text = element_text(size = 10)) + 
  facet_wrap(~Comparison, nrow =2) +
  theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(-0.04,0.2))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, label.position = "right",
                             title.position = "left", title.vjust = 0.1, legend.spacing.x = unit(0.01, 'cm'),
                            legend.spacing.y = unit(0.01, 'pt'), legend.key.size = unit(0.1, "cm"))) + 
  xlab('Signal strength')
p2 <- ggplot(pd[pd[,2]!='monocle2_trajTest',], aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=0.5)  + 
  geom_line(size=0.3) + 
  xlim(c(0,4))+
  ylab('AUC (sensitivity~realFDR)') +
  scale_color_manual(values = colv)  + 
  theme_classic() +
  theme(strip.text.x = element_text(size = 10), axis.title = element_text(size = 10), axis.text = element_text(size = 10)) + 
  facet_wrap(~Comparison, nrow =2)+
  theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(0,1))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, label.position = "right",
                             title.position = "left", title.vjust = 0.1, legend.spacing.x = unit(0.01, 'cm'),
                            legend.spacing.y = unit(0.01, 'pt'), legend.key.size = unit(0.1, "cm"))) + 
  xlab('Signal strength')
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()



