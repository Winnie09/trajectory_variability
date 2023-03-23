af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/represult')
pd <- sapply(af,function(f) {
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/represult/',f))
  d <- do.call(rbind,d)
  mean(d[,1] < 0.05) * 100
})
pd <- data.frame(sub=sub('.rds','',names(pd)),perf=pd)
library(ggplot2)
library(RColorBrewer)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/plot/perf/significance_barplot.pdf', width = 3, height = 2.5)
ggplot(pd,aes(x=sub,y=perf, fill = sub)) + 
  geom_bar(stat='identity') + 
  theme_classic()+
  theme(legend.position = 'none') +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, 'Blues')[2:9])(nrow(pd))) +
  xlab('remove cell proportion') + ylab('significance percentage')
dev.off()



