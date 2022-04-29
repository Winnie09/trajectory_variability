library(ggplot2)
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex_addBatch/'
pdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/plot/sex_addBatch/'

#rdir = '/Users/wenpinhou/Dropbox/trajectory_variability/tb/data/sex_addBatch/'
tb = read.csv(paste0(rdir, 'sample_batch_sex.csv'), row.names = 1)
pd = reshape2::melt(as.matrix(tb))
colnames(pd) = c('batch', 'sex', 'num')
pd[,1] = factor(sub('batch', '', pd[,1]),levels=sort(as.numeric(unique(sub('batch', '', pd[,1])))))

source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)

library(ggsci)
pdf(paste0(pdir, 'sample_batch_sex_barplot.pdf'), width = 5, height = 1.8)
ggplot(data = pd, aes(x = batch, y = num, fill = sex)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) + 
  scale_fill_brewer(palette = 'Pastel1') +
  ylab('number of samples') + xlab('batch')
dev.off()


