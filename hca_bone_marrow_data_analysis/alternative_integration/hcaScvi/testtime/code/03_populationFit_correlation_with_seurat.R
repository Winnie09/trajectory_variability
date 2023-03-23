library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')

df = data.frame(scvi = c('erythroid', 'myeloid', 'lymphoid'), 
                seu = c('erythroid', 'monocyte', 'lymph'))
corr <- lapply(1:3, function(j){
  print(df[j, ])
  scvi.res = readRDS(paste0('hcaScvi/testtime/res/', df[j, 1], '.rds'))
  seu.res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime/result/EM_pm/', df[j, 2], '/testtime_res.rds'))
  gene = intersect(rownames(scvi.res[[1]]), rownames(seu.res[[1]]))
  
  seu.res$populationFit <- getPopulationFit(seu.res, gene = gene, type = 'time')
  scvi.res$populationFit <- getPopulationFit(scvi.res, gene = gene, type = 'time')
  
  var1 = apply(seu.res$populationFit[gene, ], 1, var)
  var2 = apply(seu.res$populationFit[gene, ], 1, var)
  
  gene = gene[var1 > mean(var1) & var2 > mean(var2)]
  str(gene)
  
  v = sapply(gene, function(g){
    cor(seu.res$populationFit[g, ], scvi.res$populationFit[g, ])
  })
  pd = data.frame(correlation = v, lineage = df[j, 1])
  pd
})

pd = do.call(rbind, corr)
str(pd)
saveRDS(pd, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/testtime/compare/populatonFit_correlation_with_seurat_res.res')

library(ggplot2)
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/testtime/plot/populationFit_correlation_with_seurat.pdf', height = 2.1, width = 3.3)
ggplot(data = pd, aes(color = lineage, x = lineage, y = correlation)) +
  geom_density(stat = 'density') +
  labs(x="correlation of population fit", y = "density")+
  scale_color_brewer(palette = 'Dark2') + theme_classic()
dev.off()



