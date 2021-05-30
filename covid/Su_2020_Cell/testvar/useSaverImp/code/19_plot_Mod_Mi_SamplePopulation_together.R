rm(list=ls())
library(here)
# setwd(here())
setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
test.method = 'EM_pm'
comparison = 'Mod_Mi'
print(comparison)
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method, '/', comparison, '/')
pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', test.method, '/', comparison, '/')

Res <- readRDS(paste0(rdir, paste0('numeric_res_with_clu.rds')))
tb <- read.csv(paste0('/Users/wenpinhou/Dropbox/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/EM_pm/', comparison, '/differential_genes_zeyu.csv'), row.names = 1)
i = tb[1:2,1]
str(tb)

for (c in unique(tb[,11])){
  print(c)
  i <- tb[tb[,11]==c, 1]
  png(paste0(pdir, 'cluster', c, '_SamplePopulation.png'), width = 2000, height = 2000, res = 200)
  plotGeneSampleAndPopulation(Res, gene = i, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.1, continuous = F, palette = 'Dark2') 
  dev.off()  
}  

# Cluster 1—On the pre-activation stage, CD8 T cells from moderate patients are more sensitive to IFNG signaling (e.g. JAK3, IFNGR2) and pre-programed to an effector activation state by having more effector TF ZEB2 expression. This difference diminished smoothly through activation.
# Cluster 3—Along activation trajectory, CD8 T cells from moderate patients continually presents more inflammation signaling (NFKB1a, RelA) and get more inflammation driven activation (AP-1 family such as FOSB/JUNB). 
# Cluster 5—From activation to hyper-effector activation, CD8 T cells from moderate patients significantly gained more terminal effector T cell marker GNLY and secret more chemokine CCL5, with potential epigenetic changes on the global DNA methylation level (DNMT1).

glist <- list(c('TBX21', 'ZEB2'),
              c('JAK3', 'IFNGR2'),
              c('JUN', 'FOS'),
              c('CD69','IL7R'),
              c('FOSB', 'JUNB'),
              c('IKBIP','IFI6'),
              c('ISG20', 'IFT27'),
              c('IL27RA', 'FLI1'))

  pdf(paste0(pdir, 'SamplePopulation_select.pdf'), width = 3, height = 8.5)
  plotGeneSampleAndPopulation(Res, gene = unlist(glist), plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.01, line.alpha = 0.5, continuous = F, palette = 'Dark2', ncol = 2) 
  dev.off()  
  
for (c in 1:length(glist)){
  i = glist[[c]]
  png(paste0(pdir, 'cluster', c, '_SamplePopulation_select.png'), width = 600, height = 300, res = 200)
  plotGeneSampleAndPopulation(Res, gene = i, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.1, continuous = F, palette = 'Dark2') 
  dev.off()  
}

g <- c('TCF7', 'GZMB', 'IL7RA', 'KLRD1', 'PDCD1')
for (c in g){
  if (c %in% rownames(Res$covariateGroupDiff)){
    png(paste0(pdir, c, '_SamplePopulation.png'), width = 600, height = 450, res = 200)
    plotGeneSampleAndPopulation(Res, gene = c, plot.point = F, variable = 'type', axis.text.blank = T, line.size = 0.05, continuous = F, palette = 'Dark2') 
    dev.off()  
  }
    
}  





