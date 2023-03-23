rm(list=ls())
setwd('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/')
library(ggplot2)
dect.rmall <- sapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
    res <- read.csv(paste0('./res/detection.rate.csv'), row.names = 1)
  } else {
    res <- read.csv(paste0('./rmall/', pct, '/result/detection.rate.csv'), row.names = 1)
  }
    
  res['c(5, 1)', 1]
})
  
dect.rmhalf <- sapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
    res <- read.csv(paste0('./res/detection.rate.csv'), row.names = 1)
  } else {
    res <- read.csv(paste0('./rmBM1256/', pct, '/result/detection.rate.csv'), row.names = 1)
  }
  res['c(5, 1)', 1]
})
pd = rbind(data.frame(rm.pectage = seq(0, 80, 10), dect = dect.rmall, type = 'rm.all'), ## remove cells on branch 5-> 1 from all samples
            data.frame(rm.pectage = seq(0, 80, 10), dect = dect.rmhalf, type = 'rm.BM1256')) ## remove cells on branch 5-> 1 from only sample BM1256

pd[,3] <- as.factor(pd[,3])

pdf('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/perf/detectionrate.pdf', width = 2.5, height = 2.5)
ggplot(data = pd, aes(x = rm.pectage, y = dect, color = type)) +
  geom_point() +
  scale_color_brewer(palette = 'Dark2') + 
  theme_classic() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  geom_smooth(method = 'loess', level = 0, size = 0.5) +
  ylab('Detection rate') +
  xlab('Reduced cell percentage') 
dev.off()

mean.rmall <- lapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
   res <- read.csv(paste0('./res/sample.cellcomp.mean.csv'), row.names = 1)  
  } else {
   res <- read.csv(paste0('./rmall/', pct, '/result/sample.cellcomp.mean.csv'), row.names = 1)  
  }
  data.frame(samplename = paste0(colnames(res), '_rmall'), cellprop = unlist(res['c(5, 1)', ]), pct = pct*100, type = 'reduce cells on all samples', sample = 'cell reduced')
  
})
mean.rmall <- do.call(rbind, mean.rmall)
  
mean.rmhalf <- lapply(seq(0, 0.8, 0.1), function(pct){
  if (pct == 0){
   res <- read.csv(paste0('./res/sample.cellcomp.mean.csv'), row.names = 1)  
  } else {
    res <- read.csv(paste0('./rmBM1256/', pct, '/result/sample.cellcomp.mean.csv'), row.names = 1)
  }
  data.frame(samplename = paste0(colnames(res), '_rmBM1256'), cellprop = unlist(res['c(5, 1)', ]), pct = pct*100, type = 'reduce cells on half samples', sample = rep(c(rep('cell reduced',2), rep('cell not reduced',2)), 2))
})
mean.rmhalf <- do.call(rbind, mean.rmhalf) 

pdf('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/perf/cellproption.pdf', width = 4.1, height = 2.6)
ggplot(data = rbind(mean.rmall, mean.rmhalf), aes(x = pct, y = cellprop, color = sample, group = samplename)) +
  geom_line(se = F) +
  geom_point(alpha = 0.9) +
  scale_color_brewer(palette = 'Dark2') +
  facet_grid(~type) +
  theme_classic() +
  ylab('Samples\' cell proportion') +
  xlab('Reduced cell percentage')+
  theme(legend.position = 'bottom')+
  labs(color = '')
dev.off()

##-------------------------------------------------------------------
##### plot samples' proportion on myeloid lineage (5 ->1) data itself
##-------------------------------------------------------------------
cellprop.d <- readRDS('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/res/rmBM1256_cellprop_data.rds')
d = do.call(rbind, cellprop.d)
d = reshape2::melt(d) 
d.rmhalf = data.frame(samplename = paste0(d[,2], '_rmBM1256'), cellprop = d[,3], pct = as.numeric(d[,1])*100, type = 'reduce cells in half samples', sample = rep(rep(c(rep('cell reduced',2), rep('cell not reduced',2)), 2), each = length(unique(d[,1]))))

cellprop.d <- readRDS('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/res/rmall_cellprop_data.rds')
d = do.call(rbind, cellprop.d)
d = reshape2::melt(d) 
d.rmall = data.frame(samplename = paste0(d[,2], '_rmBM1256'), cellprop = d[,3], pct = as.numeric(d[,1])*100, type = 'reduce cells in all samples', sample =  'cell reduced')

pdf('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/perf/cellproption_data.pdf', width = 4.1, height = 2.6)
ggplot(data = rbind(d.rmall, d.rmhalf), aes(x = pct, y = cellprop, color = sample, group = samplename)) +
  geom_line(se = F) +
  geom_point(alpha = 0.9) +
  scale_color_brewer(palette = 'Dark2') +
  facet_grid(~type) +
  theme_classic() +
  ylab('Samples\' cell proportion') +
  xlab('Reduced cell percentage')+
  theme(legend.position = 'bottom')+
  labs(color = '')
dev.off()

