rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(grid)
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
setwd(here())
# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('function/01_function.R')
# for (path in c('monocyte','erythroid','lymph')){
for (path in c('erythroid','monocyte', 'clp', 'lymph')){
  print(path)
  ddir <- paste0('hca/real/testtime/result/EM_pm/', path)
  rdir <- pdir  <- paste0('hca/real/testtime/plot/EM_pm/', path)
  dir.create(pdir, recursive = T)
  
  # source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
  Res <- readRDS(paste0(ddir, '/testtime_res.rds'))
  ## get the differential genes
  res = Res$statistics
  res <- res[order(res[,1], -res[,3]), ]
  diffgene <- rownames(res)[res[,1] < 0.05]
  str(diffgene)
  ## --------------
  ## population fit
  ## --------------
  Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'time')
  colnames(Res$populationFit) <- colnames(Res$expr)
  saveRDS(Res$populationFit, paste0(pdir, '/population_fit.rds')) 
  
  ## -----------
  ## clustering
  ## -----------
  set.seed(12345)
  Res$cluster <- clusterGene(Res, gene = diffgene, type = 'time', k=5)
  
  ## sort clusters according to gene max expression
  v <- sapply(unique(Res$cluster), function(i){
    tmp <- Res$populationFit[Res$cluster == i, ]
    mean(apply(tmp, 1, which.max))
  })
  names(v) <- unique(Res$cluster)
  trans <- cbind(as.numeric(names(v)),rank(v))
  n <- names(Res$cluster)
  Res$cluster <- trans[match(Res$cluster,trans[,1]),2]
  names(Res$cluster) <- n
  saveRDS(Res$cluster, paste0(pdir, '/cluster.rds'))

  ## --------------
  ## save diff gene
  ## --------------
  allg <- rownames(Res$statistics[Res$statistics[,1]<0.05,,drop=FALSE])
  res <- Res$statistics[allg, ]
  res <- res[order(res[,1], -res[,3]), ]
  write.csv(res, paste0(pdir, '/testtime_differential_genes.csv'))
  
  ## ---------------
  ## plotClusterMean
  ## ----------------
  pdf(paste0(pdir, '/cluster_mean.pdf'), width = 5, height = 3.5)
  plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'time')
  dev.off()
  
  ## -----------
  ## GO analysis
  ## -----------
  goRes <- GOEnrich(testobj = Res, type = 'time', sep = ':.*')
  saveRDS(goRes, paste0(pdir, '/goRes.rds'))
  
  nn <- sapply(1:length(goRes), function(i){
    tmp <- goRes[[i]]
    # tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
    write.csv(tmp, paste0(pdir, '/cluster', i, '_GO.csv'))
    saveRDS(tmp, paste0(pdir, '/cluster', i, '_GO.rds'))
    print(str(tmp))
    return(0)
  })
  
  pdf(paste0(pdir, '/hm_GO_term.pdf'), width = 7.2, height = 3.5)
  print(plotGOEnrich(goRes))
  dev.off()
  
  # ------------------------------------------------------
  # compare original and fitted expression: not tested yet
  # ------------------------------------------------------
  
  
  png(paste0(pdir, '/fitHm.png'),width = 4000,height = 2200,res = 200)
  plotFitHm(Res, type = 'time', cellHeightTotal = 250, cellWidthTotal=200, subsampleCell=FALSE )
  dev.off()
  # png(paste0(pdir, '/fitHm_rownames.png'),width = 12000,height = 10000,res = 300)
  # print(plotFitHm(Res, showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
  # dev.off()
}
  






