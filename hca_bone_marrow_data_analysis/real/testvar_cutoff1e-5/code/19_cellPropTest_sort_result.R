rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
res <- list()
for (path in c('erythroid', 'lymph', 'monocyte')){
  print(path)
  for (variable in c('age', 'gender')){
    print(variable)
    # path = 'erythroid'
    # variable = 'age'
    pdir <- paste0('hca/real/testvar/plot/EM_pm/', path, '/', variable, '/')
    Res <- readRDS(paste0('hca/real/testvar/result/EM_pm/', path, '/', variable, '/cell_proportion_test_', variable, '.rds'))
    
    w = ifelse(variable == 'gender', 2, 1.8)
    pdf(paste0(pdir, 'cellPropTest.pdf'), width = w, height = 1.3)
    plotGene(testobj = Res, gene='prop', variable = variable,cellProp = TRUE,  variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = F, line.alpha = 1, line.size = 1, point.alpha=1, point.size=0.5, continuous = ifelse(variable=='age',T,F), sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = T, y.lab = 'Cell density')
    dev.off()
    
    res[[paste0(path, ';', variable)]] <- cbind(path=path, variable=variable, Res$statistics)    
  }
}
res <- do.call(rbind, res)
write.csv(res, 'hca/real/testvar/result/EM_pm/perf/cellprop_res.csv')

