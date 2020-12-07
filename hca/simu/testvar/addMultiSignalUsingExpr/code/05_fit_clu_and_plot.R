library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))

for (dataType  in seq(1,4)){
  ddir <- here('hca','data','simu','testvar','addMultiSignalUsingExpr')
  rdir <- here('hca','simu','testvar','addMultiSignalUsingExpr','result','addsignal', 'cluster',dataType)
  pdir <- here('hca','simu','testvar','addMultiSignalUsingExpr','plot','addsignal', dataType)
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  Res <- readRDS(here('hca','simu','testvar','addMultiSignalUsingExpr','result','addsignal', 'EM_NOT_centered',paste0(dataType,'.rds')))
  fit <- t(sapply(names(Res$fdr[Res$fdr < 0.05]), function(g){
    tmp <- get_population_fit(Res, variable = NA, gene = g)
  }))
  saveRDS(fit, paste0(rdir, '/population_fit.rds')) ##########
  
  fit.scale = scalematrix(fit) #############
  ## -----------
  ## clustering
  ## -----------
  clu.true = readRDS(paste0(ddir, '/null/geneCluster.rds'))
  set.seed(12345)
  clu <- kmeans(fit.scale, max(clu.true), iter.max = 1000)$cluster
  clu <- sort(clu)
  table(clu)
  saveRDS(clu, paste0(rdir, '/cluster.rds'))
  
  ### use population levels
  #### average the population fit for genes within a same cluster
  fitClu <- sapply(1:max(clu), function(i){
    colMeans(fit[names(clu[clu==i]), 1:(ncol(fit)),drop = FALSE])
  })
  colnames(fitClu) <- paste0('cluster', unique(clu), '(',table(clu), ')')
  
  ##### plot population fit on original scale
  pd <- reshape2::melt(fitClu)
  colnames(pd) <- c('cell', 'cluster', 'expr')
  pd$pseudotime <- Res$pseudotime[pd[,1]]
  cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
  pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
  saveRDS(pd, paste0(pdir, '/clu_based_on_population_pattern_pd.rds'))
  pdf(paste0(pdir, '/clu_based_on_population_pattern.pdf'), width = 4.5, height = 3)  
  print(ggplot(data= pd, aes(x = pseudotime, y = expr, group = cluster, color = cluster)) + 
    geom_smooth(size = 1) +
    theme_classic() +
    scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
    xlab('Pseudotime') +
    ylab('Expression') +
      theme(axis.text.x = element_blank()))
  dev.off()
  
  ##### plot population fit on standadized scale
  fitClu.scale = t(scalematrix(t(fitClu)))
  pd <- reshape2::melt(fitClu.scale)
  colnames(pd) <- c('cell', 'cluster', 'expr')
  pd$pseudotime <- Res$pseudotime[pd[,1]]
  cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
  pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
  saveRDS(pd, paste0(pdir, '/clu_based_on_population_pattern_scaled_pd.rds'))  
  pdf(paste0(pdir, '/clu_based_on_population_pattern_scaled.pdf'), width = 4.5, height = 3.2)  
  print(ggplot(data= pd, aes(x = pseudotime, y = expr, group = cluster, color = cluster)) + 
    geom_smooth(size = 1) +
    theme_classic() +
    scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
    xlab('Pseudotime') +
    ylab('Expression')) +
    theme(axis.text.x = element_blank())
  dev.off()
  
  ### use predict.values
  pred <- Res$predict.values[names(Res$fdr[Res$fdr < 0.05]), names(Res$pseudotime)]
  agg <- t(sapply(unique(clu), function(i) colMeans(pred[clu == i, ,drop=F])))
  agg <- agg[, names(Res$pseudotime)]
  rownames(agg) <- paste0('cluster', unique(clu), '(',table(clu), ')')
  
  pd <- reshape2::melt(agg)
  pd$pseudotime = Res$pseudotime[pd[,2]]
  colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
  cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
  pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
  saveRDS(pd, paste0(pdir, '/gene_cluster_pattern_pd.rds'))
  
  pdf(paste0(pdir, '/gene_cluster_pattern.pdf'), width = 4.5, height = 3.2)  
  print(ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
    geom_smooth(size = 2) +
    theme_classic() +
    scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
    xlab('Pseudotime') +
    ylab('Averaged Fitted Expression'))
  dev.off()
  
  #### scaled 
  agg.scale <- scalematrix(agg)
  pd <- reshape2::melt(agg.scale)
  pd$pseudotime = Res$pseudotime[pd[,2]]
  colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
  cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
  pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
  saveRDS(pd, paste0(pdir, '/gene_cluster_pattern_scaled_pd.rds'))
  
  pdf(paste0(pdir, '/gene_cluster_pattern_scaled.pdf'), width = 4.5, height = 3.2)  
  print(ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
    geom_smooth(size = 2) +
    theme_classic() +
    scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, 'Set1')))(max(clu))) +
    xlab('Pseudotime') +
    ylab('Scaled Averaged Fitted Expression'))
  dev.off()
  
  ## ----------------
  ## confusion matrix
  ## ----------------
  library(caret)
  fromgene = readRDS(paste0(ddir, '/fromgene/', dataType, '.rds'))
  selgene = readRDS(paste0(ddir, '/selgene/selgene.rds'))
  clu.true = clu.true[fromgene]
  clu.true = paste0('cluster', clu.true)
  names(clu.true) <- selgene
  
  m <- confusionMatrix(data = as.factor(paste0('cluster',clu)), reference = as.factor(clu.true[names(clu)]), 
                       dnn = c('diffGeneCluster','SignalCluster'))  ## dnn = c("Prediction", "Reference")
  tb <- t(m$table)
  tb <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm.pdf'), width = 3.5, height = 3)
  print(pheatmap(tb))
  dev.off()
}
  

