library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')

df = data.frame(scvi = c('erythroid', 'myeloid', 'lymphoid'), 
                seu = c('erythroid', 'monocyte', 'lymph'))
corr <- lapply(1:3, function(j){
  print(df[j, ])
  scvi.res = readRDS(paste0('hcaScvi/testvar/res/', df[j, 1], '.rds'))
  seu.res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/', df[j, 2], '/gender/gender_res.rds'))
  gene = intersect(rownames(scvi.res[[1]]), rownames(seu.res[[1]]))
  str(gene)
  seu.res$populationFit <- getPopulationFit(seu.res, gene = gene, type = 'variable')
  
  colnames(scvi.res$design) = c('intercept', 'gender') 
  scvi.res$populationFit <- getPopulationFit(scvi.res, gene = gene, type = 'variable')
  
  var1 = apply(seu.res$populationFit[[1]][gene, ], 1, var)
  var2 = apply(seu.res$populationFit[[2]][gene, ], 1, var)
  gene = gene[var1 > mean(var1) & var2 > mean(var2)]
  str(gene)
  
  v = t(sapply(gene, function(g){
    c(cor(seu.res$populationFit[[1]][g, ], scvi.res$populationFit[[2]][g, ]),
      cor(seu.res$populationFit[[2]][g, ], scvi.res$populationFit[[1]][g, ]))  
  }))
  pd = reshape2::melt(data.frame(female = v[,1], male = v[,2], lineage = df[j, 1]))
  colnames(pd) = c('lineage', 'sex', 'correlation')
  pd
})

pd = do.call(rbind, corr)
str(pd)
saveRDS(pd, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/testvar/compare/populatonFit_correlation_with_seurat_res.res')

library(ggplot2)
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/testvar/plot/populationFit_correlation_with_seurat.pdf', height = 2.4, width = 3.6)
ggplot(data = pd, aes(color = lineage, x = correlation)) +
  geom_density(aes(linetype = sex)) +
  labs(x="correlation of population fit", y = "density")+
  scale_color_brewer(palette = 'Dark2') + theme_classic()
dev.off()



