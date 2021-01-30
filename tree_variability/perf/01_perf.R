rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/')
library(ggplot2)
dect.rmall <- sapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
    res <- read.csv(paste0('./auto_pc_auto_nclu_module/result/detection.rate.csv'), row.names = 1)
  } else {
    res <- read.csv(paste0('./rmall/', pct, '/result/detection.rate.csv'), row.names = 1)
  }
    
  res['c(5, 1)', 1]
})
  
dect.rmhalf <- sapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
    res <- read.csv(paste0('./auto_pc_auto_nclu_module/result/detection.rate.csv'), row.names = 1)
  } else {
    res <- read.csv(paste0('./rmBM1256/', pct, '/result/detection.rate.csv'), row.names = 1)
  }
  res['c(5, 1)', 1]
})
pd = rbind(data.frame(rm.pectage = seq(0, 80, 10), dect = dect.rmall, type = 'rm.all'),
            data.frame(rm.pectage = seq(0, 80, 10), dect = dect.rmhalf, type = 'rm.BM1256'))

pd[,3] <- as.factor(pd[,3])

pdf('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/detectionrate.pdf', width = 2.3, height = 2.3)
ggplot(data = pd, aes(x = rm.pectage, y = dect, color = type)) +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') + 
  theme_classic() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  geom_smooth() +
  ylab('Detection Rate') +
  xlab('Removed Percentage') 
dev.off()

mean.rmall <- lapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
   res <- read.csv(paste0('./auto_pc_auto_nclu_module/result/sample.cellcomp.mean.csv'), row.names = 1)  
  } else {
   res <- read.csv(paste0('./rmall/', pct, '/result/sample.cellcomp.mean.csv'), row.names = 1)  
  }
  data.frame(sample = paste0(colnames(res), '_rmall'), cellprop = unlist(res['c(5, 1)', ]), pct = pct*100, type = 'rm.all', sample.cells = rep(c(rep('Removed',2), rep('Unremoved',2)), 2))
  
})
mean.rmall <- do.call(rbind, mean.rmall)
  
mean.rmhalf <- lapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
   res <- read.csv(paste0('./auto_pc_auto_nclu_module/result/sample.cellcomp.mean.csv'), row.names = 1)  
  } else {
    res <- read.csv(paste0('./rmBM1256/', pct, '/result/sample.cellcomp.mean.csv'), row.names = 1)
  }
  data.frame(sample = paste0(colnames(res), '_rmBM1256'), cellprop = unlist(res['c(5, 1)', ]), pct = pct*100, type = 'rm.BM1256', sample.cells = rep(c(rep('Removed',2), rep('Unremoved',2)), 2))
})
mean.rmhalf <- do.call(rbind, mean.rmhalf) 

pdf('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/cellproption.pdf', width = 5, height = 2.2)
ggplot(data = rbind(mean.rmall, mean.rmhalf)) +
  geom_smooth(aes(x = pct, y = cellprop, color = sample.cells, group = sample), se = F) +
  scale_color_brewer(palette = 'Dark2') +
  facet_grid(~type) +
  theme_classic() +
  ylab('Cell Proportion') +
  xlab('Removed Cells Percentage')
dev.off()

