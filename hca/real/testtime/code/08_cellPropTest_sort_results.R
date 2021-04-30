library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
res <- list()
for (path in c('erythroid', 'lymph', 'monocyte')){
    pdir <- paste0('hca/real/testtime/plot/EM_pm/', path, '/')
    Res <- readRDS(paste0('hca/real/testtime/result/EM_pm/', path, '/cell_proportion_test.rds'))
    
    pdf(paste0(pdir, 'cellPropTest.pdf'), width = 3, height = 2)
    plotGene(testobj = Res, gene='prop', cellProp = TRUE,  variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = F, line.alpha = 1, line.size = 1, point.alpha=1, point.size=0.5, continuous = T, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = T, y.lab = 'Cell density')
    dev.off()
    
    pt = Res$pseudotime
    expr = Res$expr 
    cellanno = Res$cellanno
    pred <- predict_fitting(Res, gene = 'prop', test.type = Res$test.type)
    ld <- data.frame(den = pred[1,], pt = pt[colnames(pred)], sample = cellanno[match(colnames(pred),cellanno[,1]),2])
    pd <- data.frame(den = Res$expr[1,], pt = pt[colnames(Res$expr)], sample = cellanno[match(colnames(expr), cellanno[,1]), 2])
    library(ggplot2)
    pdf(paste0(pdir, 'cellPropTest_fitting.pdf'), width = 3, height = 6)
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
write.csv(res, 'hca/real/testtime/result/EM_pm/perf/cellprop_res.csv')

