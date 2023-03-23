ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/cellprop/result/testres/'
af <- list.files(ddir)
pd <- sapply(af,function(f) {
  print(f)
  d <- readRDS(paste0(ddir, f))
  d <- do.call(rbind,d)
  mean(d[complete.cases(d),1] < 0.05) * 100
})
pd <- data.frame(sub=sub('.rds','',names(pd)),perf=pd)
library(ggplot2)
library(RColorBrewer)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/cellprop/plot/perf/significance_barplot.pdf', width = 3.3, height = 2.5)
ggplot(pd,aes(x=sub,y=perf, fill = sub)) + 
  geom_bar(stat='identity') + 
  theme_classic()+
  theme(legend.position = 'none') +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, 'Blues')[3:9])(nrow(pd))) +
  xlab('remove cell proportion') + ylab('significance percentage')
dev.off()

