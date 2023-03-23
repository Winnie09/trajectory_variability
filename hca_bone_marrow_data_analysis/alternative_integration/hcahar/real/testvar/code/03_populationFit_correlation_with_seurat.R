library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
lin = c('erythroid', 'monocyte', 'lymph')
corr <- lapply(lin, function(j){
  har.res = readRDS(paste0('hcahar/real/testvar/result/', j, '/EM_pm.rds'))
  seu.res = readRDS(paste0('hca/real/testvar/result/EM_pm/', j, '/gender/gender_res.rds'))
  
  gene = intersect(rownames(har.res[[1]]), rownames(seu.res[[1]]))
  str(gene)
  seu.res$populationFit <- getPopulationFit(seu.res, gene = gene, type = 'variable')
  
  colnames(har.res$design) = c('intercept', 'gender') 
  har.res$populationFit <- getPopulationFit(har.res, gene = gene, type = 'variable')
  
  var1 = apply(seu.res$populationFit[[1]][gene, ], 1, var)
  var2 = apply(seu.res$populationFit[[2]][gene, ], 1, var)
  gene = gene[var1 > mean(var1) & var2 > mean(var2)]
  str(gene)
  
  v = t(sapply(gene, function(g){
    c(cor(seu.res$populationFit[[1]][g, ], har.res$populationFit[[2]][g, ]),
      cor(seu.res$populationFit[[2]][g, ], har.res$populationFit[[1]][g, ]))  
  }))
  pd = reshape2::melt(data.frame(female = v[,1], male = v[,2], lineage = j))
  colnames(pd) = c('lineage', 'sex', 'correlation')
  pd
})

pd = do.call(rbind, corr)
str(pd)
dir.create('hcahar/testvar/compare')
saveRDS(pd, 'hcahar/testvar/compare/populatonFit_correlation_with_seurat_res.res')

library(ggplot2)
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('hcahar/real/testvar/plot/populationFit_correlation_with_seurat.pdf', height = 2.4, width = 3.6)
ggplot(data = pd, aes(color = lineage, x = correlation)) +
  geom_density(aes(linetype = sex)) +
  labs(x="correlation of population fit", y = "density")+
  scale_color_brewer(palette = 'Dark2') + theme_classic()
dev.off()




