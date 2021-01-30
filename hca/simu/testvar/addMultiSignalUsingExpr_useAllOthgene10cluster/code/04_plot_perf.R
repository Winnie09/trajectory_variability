rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/perf/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
pd <- readRDS(paste0(rdir, 'perf.rds'))

pd[,2] = as.factor(pd[,2])
pd[,1] = as.numeric(pd[,1])
pd[,3] = as.numeric(pd[,3])
pd[,4] = as.numeric(pd[,4])
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
pdf(paste0(pdir, 'fdr_diff_auc.pdf'), width =11, height = 3.2)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=1.5)  + 
  geom_line(size=0.2) + 
  theme_classic() + 
  ylab('fdr.diff(real~reported - 0.25*0.25/2)') +
  scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(length(unique(pd$Method))))
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=1.5)  + 
  geom_line(size=0.2) + 
  theme_classic() + 
  ylab('Area Under Sensitivity Real_FDR curve') +
  scale_color_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(length(unique(pd$Method))))
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()

