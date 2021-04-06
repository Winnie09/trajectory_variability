rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
ddir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
# for (comparison in c('HD_Se', 'HD_Mi','Mod_Se','HD_Mod', 'Se_Mi', 'Recovered_Deceased')){
for (comparison in c('Mod_Mi', 'Se_Mi', 'Recovered_Deceased')){
  print(comparison)
  rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', comparison, '/')
  pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', comparison, '/')
  dir.create(rdir, recursive = T)
  dir.create(pdir, recursive = T)
  Res <- readRDS(paste0(ddir, paste0('numeric_', comparison, '_res.rds')))
  print(names(Res))


  diffgene <- rownames(Res$statistics[Res$statistics[,7]<0.05, ])
  str(diffgene)
  
  ## --------------
  ## population fit
  ## --------------
  Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')
  
  ## -----------
  ## clustering
  ## -----------
  Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene)
  Res$cluster <- clusterGene(Res, gene = diffgene, type = 'variable', k=3)
  
  ## --------------
  ## save diff gene
  ## --------------
  DEGType <- getDEGType(Res)
  Res$DEGType <- DEGType
  res <- Res$statistics
  res <- data.frame(res, DEGType = DEGType[rownames(res)], stringsAsFactors = F)
  saveRDS(res, paste0(rdir, 'differential_genes.rds'))
  write.csv(res, paste0(rdir, 'differential_genes.csv'))

  ## ----------------
  ## plotClusterDiff
  ## ----------------
  pdf(paste0(pdir, '/cluster_diff.pdf'), width = 3, height = 2)
  plotClusterDiff(testobj=Res, gene = diffgene)
  dev.off()
  
  ## ---------------
  ## plotClusterMean
  ## ----------------
  pdf(paste0(pdir, '/cluster_mean.pdf'), width = 5, height = 3.5)
  plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'variable')
  dev.off()
    
  ## -----------
  ## GO analysis
  ## -----------
  goRes <- GOEnrich(testobj = Res, type = 'variable')
  nn <- sapply(1:length(goRes), function(i){
    tmp <- goRes[[i]]
    tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
    if (nrow(tmp)>0) {write.csv(tmp, paste0(rdir, 'cluster', i, '_GO.csv'))}
    print(str(tmp))
    return(0)
  })
  
  pdf(paste0(pdir, 'hm_GO_term.pdf'), width = 7.2, height = 3.5)
  print(plotGOEnrich(goRes))
  dev.off()

  ## -----------------------
  ## plotClusterMeanAndDiff
  ## -----------------------
  pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 4, height = 5.5)
  print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
  dev.off()
  
  res <- res[names(Res$cluster), ]
  res <- cbind(res, cluster = Res$cluster)
  for (i in 1:max(Res$cluster)){
    print(i)
    gene <- rownames(res)[res$cluster == i]
    png(paste0(pdir, 'diffgene_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
    print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
    dev.off()
  }

  for (i in 1:max(Res$cluster)) {
    print(i)
    gene <- rownames(res)[res$cluster == i]
    png(paste0(pdir, 'diffgene_groupDiff_cluster', i, '.png'),width = 2500,height = 2500, res = 200)
    print(plotClusterDiff(testobj = Res, gene = gene[1:min(length(gene), 100)], each = TRUE, sep = ':.*'))
    dev.off()
  }

  # --------------------------------------
  # compare original and fitted expression
  # --------------------------------------
  png(paste0(pdir, 'fitHm.png'),width = 4500, height = 2500,res = 300)
  print(plotFitHm(Res, type = 'variable',numSubsampleCell=1e3))
  dev.off()

  png(paste0(pdir, 'fitHm_rownames.png'),width = 11080,height = length(Res$cluster)*25,res = 200)
  print(plotFitHm(Res, showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 8, type = 'variable'))
  dev.off()
}
  

