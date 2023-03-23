rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
ddir <- 'hca/simu/testvar/nullsimu/result/perf/'
pdir <- 'hca/simu/testvar/nullsimu/plot/perf/'
dir.create(pdir, recursive = TRUE)
num <- readRDS(paste0(ddir, 'number_of_false_positives.rds'))
num = num[names(num)!='Limma']
names(num)[2] <- 'ourmethod'
pd <- data.frame(num = num, method = names(num), stringsAsFactors = FALSE)
library(ggplot2)
library(RColorBrewer)
pdf(paste0(pdir, 'num_false_positives.pdf'), height = 2.6, width=1.2)
ggplot(data = pd, aes(y = num, x = method, fill = method)) + 
  geom_bar(stat = 'identity') +
  theme_classic() +
  scale_fill_brewer(palette = 'Dark2')+
  ylab('number of false positives') +
  xlab('method')+
  theme(legend.position = 'none', 
        axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = c('red', 'black'))) 
  
dev.off()


