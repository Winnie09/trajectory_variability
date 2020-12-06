pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/plot/addsignal/'
perf <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/addsignal/perf/perf.rds')
#perf <- perf[perf[,2]!='EM_centered',]
pd$SignalStreghth <- as.numeric(pd$Type)
pd[,2] = as.factor(pd[,2])
pd[,5] = as.factor(pd[,5])
pd[,3] = as.numeric(pd[,3])
pd[,4] = as.numeric(pd[,4])
library(ggplot2)
library(gridExtra)
pdf(paste0(pdir, 'compare_fdr_diff_auc.pdf'),width=10,height=3.3)
p1 <- ggplot(pd, aes(x = SignalStreghth, y = Fdr.Diff, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('fdr.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Dark2')  
p2 <- ggplot(pd, aes(x = SignalStreghth, y = AUC, color=Method)) + 
  geom_point(size=3)  + 
  geom_line(size=1) + 
  theme_classic() + 
  ylab('Area Under Sensitivity-RealFDR curve') +
  scale_color_brewer(palette = 'Dark2') 
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

