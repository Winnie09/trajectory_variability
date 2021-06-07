af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/represult')
pd <- sapply(af,function(f) {
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/represult/',f))
  d <- do.call(rbind,d)
  d[,1]
})
colnames(pd) <- sub('.rds','',colnames(pd))
df = reshape2::melt(pd)
df[,2] = as.factor(df[,2])
library(ggplot2)
library(RColorBrewer)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/plot/perf/pvalue_violinplot.pdf', width = 3, height = 2.5)
ggplot(df, aes(x=Var2,y=value, color = Var2)) + 
  geom_violin(scale = 'width') +
  geom_jitter(size = 0.01) +
  theme_classic()+
  theme(legend.position = 'none') +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, 'Blues')[2:9])(ncol(pd))) +
  xlab('remove cell proportion') + ylab('significance percentage')
dev.off()

