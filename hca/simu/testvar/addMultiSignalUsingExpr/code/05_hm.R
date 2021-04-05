rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# source('/Users/wenpinhou/Dropbox/resource/myfunc/01_function.R')
source('function/01_function.R')
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr'
selgene <- readRDS(paste0(ddir, '/data/selgene/selgene.rds'))
selgene1 <- readRDS(paste0(ddir, '/data/selgene/selgene1.rds'))
selgene2 <- readRDS(paste0(ddir, '/data/selgene/selgene2.rds'))
selgene3 <- readRDS(paste0(ddir, '/data/selgene/selgene3.rds'))

for (dataType in seq(4, 1)) {
  print(dataType)
  pdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/plot/EM_pm/', dataType)

  dir.create(pdir, showWarnings = F, recursive = T)
  Res <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/EM_pm/',dataType,'.rds'))
  names(Res)
  res <- Res$statistics
  res <- res[order(res[,1], -res[,3]), ]
  write.csv(res, paste0(pdir, '/Res_statistics.csv'))  
  # clu <- readRDS(paste0(rdir, '/cluster.rds'))
  
  diffgene = rownames(res[res[,1] < 0.05, ])
  fit <- getPopulationFit(testobj = Res, gene = diffgene, type = Res$test.type)
  saveRDS(fit, paste0(pdir, '/population_fit.rds')) ##########
  
  # fit <- readRDS(paste0(rdir, '/population_fit.rds'))
  # str(fit)
  Res$populationFit <- fit
  fit.bak = fit
  Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene)
  
  Res$cluster <- clusterGene(Res, gene = diffgene, type = 'variable', k=7)
  clu <- Res$cluster
  saveRDS(clu, paste0(pdir, '/cluster.rds')) ##########
  saveRDS(Res, paste0(pdir, '/Res_with_clu.rds'))
  
  DDGType <- getDDGType(Res)
  DDGType <- DDGType[names(clu)]
  fit.scale <- do.call(cbind, fit)
  fit.scale <- fit.scale[names(clu), ]
  fit.scale <- scalematrix(fit.scale)
  str(fit.scale)
  colnames(fit.scale) <- seq(1, ncol(fit.scale))
  
  res <- data.frame(clu = clu, 
                    DDGType = DDGType[names(clu)],
                    cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(1, ncol(fit.scale)/2)], seq(1, ncol(fit.scale)/2))))
  res$signalType <- sapply(rownames(res), function(i) {
    if (i %in% selgene1){
      '1'
    } else if (i %in% selgene2){
      '2'
    } else {
      '3'
    }
  })
  write.csv(res, paste0(pdir, '/differential_genes.csv'))
  fit.scale <- fit.scale[rownames(res)[order(res$clu, res$signalType, res$cor)], ]
  # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
  meanres <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/meandiff/', dataType, '.rds'))
  str(meanres)
  meanDiffTest <- ifelse(meanres[rownames(fit.scale), 5] < 0.05, 'Diff', 'nonDiff')
  names(meanDiffTest) <- rownames(fit.scale)
  
  rowann = data.frame(
    cluster = as.character(clu),
    signalType = sapply(names(clu), function(i){
      if (i %in% selgene1){
        'trend only'
      } else if (i %in% selgene2){
        'mean only'
      } else {
        'both'
      }
    }),
    DDGType = DDGType[names(clu)], 
    limmaPb = meanDiffTest[names(clu)],
    gs = sapply(names(clu), function(i)
      ifelse(i %in% selgene, 'Yes', 'No')),
    stringsAsFactors = F
  )
  rownames(rowann) = names(clu)
  rowann <- rowann[rownames(fit.scale), ]
  saveRDS(rowann, paste0(pdir, '/rowann.rds'))
  
  png(
    paste0(pdir, '/hm_kmeans_population_difffit_scale.png'),
    bg = 'white',
    width = 4500,
    height = 3500,
    res = 300
  )
  plotDiffFitHm(Res, type='variable', rowann = rowann)
  dev.off()
}


# ## plot trendonly signal & limman cannot identifieid genes
# png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_genes.png'), width = 3000, height = 3000, res = 200)
# plotGene(Res, id[1:100], plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
# dev.off()
# png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_genes_population.png'), width = 3000, height = 3000, res = 200)
# plotGenePopulation(Res, id[1:100], type = 'variable', sep = ':.*')
# dev.off()
# 
# ## plot selected example genes
# allg <- c('CTDP1', 'CFDP1', 'NOL6', 'TMEM198B')
# allg <- sapply(allg, function(i){
#   rownames(fit.scale)[grepl(i, rownames(fit.scale))]
# })
# png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_example.png'), width = 800, height = 550, res = 200)
# plotGene(Res, allg, plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
# dev.off()
# png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_example_population.png'), width = 850, height = 500, res = 200)
# plotGenePopulation(Res, allg, type = 'variable', sep = ':.*')
# dev.off()
# 
# ## plot cluster mean and difference
# pdf(paste0(pdir, '/clusterMeanAndDiff.pdf'), width = 3, height = 6)
# plotClusterMeanAndDiff(Res, clu)
# dev.off()
# 
# ## plot diffType == unknown gene
# png(paste0(pdir, '/diffType_unknown_gene.png'), width = 5000, height = 5000, res = 200)
# plotGene(Res, names(diffType[diffType=='unknown'])[1:100], plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
# dev.off()
# 
# png(paste0(pdir, '/diffType_meanOnly_gene.png'), width = 5000, height = 5000, res = 200)
# plotGene(Res, names(diffType[diffType=='meanOnly'])[1:100], plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
# dev.off()
# 
# png(paste0(pdir, '/diffType_trendOnly_gene.png'), width = 5000, height = 5000, res = 200)
# plotGene(Res, names(diffType[diffType=='trendOnly'])[1:100], plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
# dev.off()
# 
# png(paste0(pdir, '/diffType_both_gene.png'), width = 5000, height = 5000, res = 200)
# plotGene(Res, names(diffType[diffType=='both'])[1:100], plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
# dev.off()
# 
# 
# 
# 
# 

