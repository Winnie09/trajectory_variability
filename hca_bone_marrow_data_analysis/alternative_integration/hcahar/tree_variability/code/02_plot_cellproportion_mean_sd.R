rm(list=ls())
library(gridExtra)
library(RColorBrewer)
setwd('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/')
library(ggplot2)
res <- read.csv(paste0('./res/sample.cellcomp.mean.csv'), row.names = 1)
pd1 <- reshape2::melt(as.matrix(res))
colnames(pd1) <- c('Branch', 'Sample', 'Cell.proportion')
p1 <-ggplot(data = pd1, aes(x = Branch, y = Sample, fill = Cell.proportion)) + geom_tile() + 
  scale_fill_distiller('Cell proportion (mean)',palette = "RdYlBu") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

res <- read.csv(paste0('./res/sample.cellcomp.sd.csv'), row.names = 1)
pd2 <- reshape2::melt(as.matrix(res))
colnames(pd2) <- c('Branch', 'Sample', 'Cell.proportion')
p2 <- ggplot(data = pd2, aes(x = Branch, y = Sample, fill = Cell.proportion)) + geom_tile() + 
  scale_fill_distiller('Cell proportion (sd)',palette = "RdYlBu", trans = 'reverse') + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



res <- read.csv(paste0('./rmall/0.8/result/sample.cellcomp.mean.csv'), row.names = 1)
pd3 <- reshape2::melt(as.matrix(res))
colnames(pd3) <- c('Branch', 'Sample', 'Cell.proportion')
p3 <- ggplot(data = pd3, aes(x = Branch, y = Sample, fill = Cell.proportion)) + geom_tile() + 
  scale_fill_distiller('Cell proportion (mean)',palette = "RdYlBu") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

res <- read.csv(paste0('./rmall/0.8/result/sample.cellcomp.sd.csv'), row.names = 1)
pd4 <- reshape2::melt(as.matrix(res))
colnames(pd4) <- c('Branch', 'Sample', 'Cell.proportion')
p4 <- ggplot(data = pd4, aes(x = Branch, y = Sample, fill = Cell.proportion)) + geom_tile() + 
  scale_fill_distiller('Cell proportion (sd)',palette = "RdYlBu", trans = 'reverse') + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/perf/cellproportion_mean_sd.pdf', width = 7, height= 4.5)
grid.arrange(p1,p2,p3,p4, nrow = 2)
dev.off()


