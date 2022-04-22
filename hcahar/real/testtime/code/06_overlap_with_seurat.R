library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
lin = c('erythroid', 'monocyte', 'lymph')
ov <- sapply(lin, function(j){
  har.res = readRDS(paste0('hcahar/real/testtime/result/', j, '/EM_pm.rds'))
  seu.res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime/result/EM_pm/', j, '/testtime_res.rds'))
  har.s = har.res[[1]]
  seu.s = seu.res[[1]]
  har.s = har.s[order(har.s[1], -har.s[3]), ]
  seu.s = seu.s[order(seu.s[1], -seu.s[3]), ]
  tmp = sapply(1:50, function(i){
    o = intersect(rownames(seu.s)[1:(i*100)], rownames(har.s)[1:(i*100)])
    length(o)/(i*100)
  })  
})
colnames(ov) = lin
rownames(ov) = paste0('top', (1:50)*100)
saveRDS(ov, 'hcahar/real/testtime/compare/overlap_with_seurat_res.res')

pd = reshape2::melt(ov)
pd[,1] = as.numeric(sapply(pd[,1], function(i) gsub('top', '', i))) 
colnames(pd) = c('num', 'lineage', 'prop')
 
library(ggplot2)
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('hcahar/real/testtime/plot/overlap_with_seurat.pdf', height = 1.4, width = 2.1)
ggplot(data = pd, aes(x = num, y = prop, color = lineage)) + 
  geom_point(size = 0.3) + 
  geom_line(size = 0.2, alpha = 0.5) + scale_color_brewer(palette = 'Dark2') +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5)) + 
  xlab('number of top genes') + 
  ylab('overlap proportion') 
dev.off()

