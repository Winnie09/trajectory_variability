library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
lin = c('erythroid', 'monocyte', 'lymph')
corr <- lapply(lin, function(j){
  seu.res = readRDS(paste0('hca/real/testtime/result/EM_pm/', j, '/testtime_res.rds'))
  har.res = readRDS(paste0('hcahar/real/testtime/result/', j, '/EM_pm.rds'))
  gene = intersect(rownames(har.res[[1]]), rownames(seu.res[[1]]))
  seu.res$populationFit <- getPopulationFit(seu.res, gene = gene, type = 'time')
  har.res$populationFit <- getPopulationFit(har.res, gene = gene, type = 'time')
  
  var1 = apply(seu.res$populationFit[gene, ], 1, var)
  var2 = apply(seu.res$populationFit[gene, ], 1, var)
  summary(var1)
  summary(var2)
  
  gene = gene[var1 > mean(var1) & var2 > mean(var2)]
  str(gene)
  
  v = sapply(gene, function(g){
    cor(seu.res$populationFit[g, ], har.res$populationFit[g, ])
  })
  pd = data.frame(correlation = v, lineage = j)
  pd
})

pd = do.call(rbind, corr)
str(pd)
saveRDS(pd, 'hcahar/real/testtime/compare/populatonFit_correlation_with_seurat_res.res')

library(ggplot2)
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('hcahar/real/testtime/plot/populationFit_correlation_with_seurat.pdf', height = 2.1, width = 3.3)
ggplot(data = pd, aes(color = lineage,  x = correlation)) +
  geom_density(stat = 'density', adjust = 3) +
  labs(x="correlation of population fit", y = "density")+
  scale_color_brewer(palette = 'Dark2') + theme_classic()
dev.off()




