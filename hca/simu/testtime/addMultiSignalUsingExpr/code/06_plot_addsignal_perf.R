setwd("/Users/wenpinhou/Dropbox/trajectory_variability")
# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
pdir <- 'hca/simu/testtime/addMultiSignalUsingExpr/plot/addsignal/perf/'
dir.create(pdir, recursive=T, showWarnings=F)
pd <- readRDS('hca/simu/testtime/addMultiSignalUsingExpr/result/addsignal/perf/perf.rds')
library(ggplot2)
library(gridExtra)
pdf(paste0(pdir, 'compare_fdr_diff_auc_inc_chisq.pdf'),width=10,height=3.2)
p1 <- ggplot(pd, aes(x = SignalStreghth, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('fdr.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Set1')  
p2 <- ggplot(pd, aes(x = SignalStreghth, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('Area Under Sensitivity-RealFDR curve') +
  scale_color_brewer(palette = 'Set1') 
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

pd <- pd[pd$Method != 'EM_chisq', ]
pdf(paste0(pdir, 'compare_fdr_diff_auc.pdf'),width=10,height=3.2)
p1 <- ggplot(pd, aes(x = SignalStreghth, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('fdr.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Set1')  
p2 <- ggplot(pd, aes(x = SignalStreghth, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('Area Under Sensitivity-RealFDR curve') +
  scale_color_brewer(palette = 'Set1') 
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()



