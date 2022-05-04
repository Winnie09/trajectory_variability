#setwd('/Users/wenpinhou/Dropbox/')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
pd = readRDS(file='trajectory_variability/tb/res/sexprop/pc2_pd.rds')
source('resource/ggplot_theme.R')
theme_set(.new_theme)
library(ggplot2)
pdf('trajectory_variability/tb/res/sexprop/pc2.pdf', width = 2.5, height = 1.6)
ggplot(data = pd,aes(pt,group=sample,col=sex)) + 
  geom_density(size = 0.1) + 
  scale_color_manual(values=c('orange','royalblue')) + 
  xlab('pseudotime') +
  ylab('cell density')
dev.off()



