library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
res <- list()
for (path in c('erythroid', 'lymph', 'monocyte')){
  print(path)
    pdir <- paste0('hca/real/testtime/plot/EM_pm/', path, '/')
    Res <- readRDS(paste0('hca/real/testtime/result/EM_pm/', path, '/gender/cell_proportion_test.rds'))
    
    pdf(paste0(pdir, 'cellPropTest_gender.pdf'), width = 2.6, height = 1.9)
    plotGene(testobj = Res, gene='prop', cellProp = TRUE,  variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = F, line.alpha = 1, line.size = 0.4, point.alpha=0.8, point.size=0.1, continuous = T, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = T, y.lab = 'Cell density')
    dev.off()
    
    pdf(paste0(pdir, 'cellPropTest_gender_point.pdf'), width = 2.6, height = 1.9)
    plotGene(testobj = Res, gene='prop', cellProp = TRUE,  variable = 'gender', variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = T, line.alpha = 1, line.size = 0.4, point.alpha=0.8, point.size=0.5, continuous = F, sep = NA, palette = 'Dark2', use.palette=T, ncol = NULL,  axis.text.blank = F, y.lab = 'Cell density')
    dev.off()
    
    pt = Res$pseudotime
    expr = Res$expr 
    cellanno = Res$cellanno
    pred <- predict_fitting(Res, gene = 'prop', test.type = Res$test.type)
    ld <- data.frame(den = pred[1,], pt = pt[colnames(pred)], sample = cellanno[match(colnames(pred),cellanno[,1]),2])
    pd <- data.frame(den = Res$expr[1,], pt = pt[colnames(Res$expr)], sample = cellanno[match(colnames(expr), cellanno[,1]), 2])
    library(ggplot2)
    pdf(paste0(pdir, 'cellPropTest_fitting_gender.pdf'), width = 3, height = 6)
    print(ggplot() + 
      geom_line(data = ld, aes(x = pt, y = den), size = 0.5)  +
      geom_point(data = pd, aes(x =pt, y = den), size = 0.2)+
      facet_wrap(~sample, nrow = 4) +
      theme_minimal() + xlab('Pseudotime bin') + ylab('Density'))
    dev.off()
    res[[path]] <- cbind(path=path, Res$statistics)    
}
res <- do.call(rbind, res)
dir.create('hca/real/testtime/result/EM_pm/perf/', recursive = T)

res[,'fdr'] = p.adjust(res[,'pval.overall'], 'fdr')
write.csv(res, 'hca/real/testtime/result/EM_pm/perf/cellprop_res.csv')
saveRDS(res, 'hca/real/testtime/result/EM_pm/perf/cellprop_res.rds')

###############################
for (path in c('erythroid', 'lymph', 'monocyte')){
  print(path)
    pdir <- paste0('hca/real/testtime/plot/EM_pm/', path, '/')
    Res <- readRDS(paste0('hca/real/testtime/result/EM_pm/', path, '/age/cell_proportion_test.rds'))
    
    pdf(paste0(pdir, 'cellPropTest_age.pdf'), width = 2.8, height = 1.9)
    plotGene(testobj = Res, gene='prop', cellProp = TRUE,  variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = T, line.alpha = 1, line.size = 0.4, point.alpha=0.8, point.size=0.5, continuous = T, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = T, y.lab = 'Cell density')
    dev.off()
    
    pdf(paste0(pdir, 'cellPropTest_age_point.pdf'), width = 2.3, height = 1.9)
    plotGene(testobj = Res, gene='prop', cellProp = TRUE,  variable = 'age', variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = T, line.alpha = 1, line.size = 0.4, point.alpha=0.8, point.size=0.5, continuous = T, sep = NA, palette = 'Dark2', use.palette=F, ncol = NULL,  axis.text.blank = F, y.lab = 'Cell density')
    dev.off()
    
    pt = Res$pseudotime
    expr = Res$expr 
    cellanno = Res$cellanno
    pred <- predict_fitting(Res, gene = 'prop', test.type = Res$test.type)
    ld <- data.frame(den = pred[1,], pt = pt[colnames(pred)], sample = cellanno[match(colnames(pred),cellanno[,1]),2])
    pd <- data.frame(den = Res$expr[1,], pt = pt[colnames(Res$expr)], sample = cellanno[match(colnames(expr), cellanno[,1]), 2])
    library(ggplot2)
    pdf(paste0(pdir, 'cellPropTest_fitting_age.pdf'), width = 3, height = 6)
    print(ggplot() + 
      geom_line(data = ld, aes(x = pt, y = den), size = 0.5)  +
      geom_point(data = pd, aes(x =pt, y = den), size = 0.2)+
      facet_wrap(~sample, nrow = 4) +
      theme_minimal() + xlab('Pseudotime bin') + ylab('Density'))
    dev.off()
}

