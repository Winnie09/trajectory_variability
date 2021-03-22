rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
ddir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
test.method = 'EM_pm'
# for (comparison in c('HD_Se', 'HD_Mi','Mod_Se','HD_Mod', 'Se_Mi', 'Recovered_Deceased')){
# for (comparison in list.files(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method))){
for (comparison in c('Se_Mi','Mod_Mi',  'Recovered_Deceased')){
  print(comparison)
  rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method, '/', comparison, '/')
  pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', test.method, '/', comparison, '/')
  
  dir.create(pdir, recursive = T, showWarnings = F)
  Res <- readRDS(paste0(rdir, paste0('numeric_res.rds')))
  print(names(Res))
  
  statistics = Res$statistics
  diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
  str(diffgene)
  ## if the above test works well, then refine the following codes
  ## --------------
  ## population fit
  ## --------------
  
  Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')
  
  ## -----------
  ## clustering
  ## -----------
  if (comparison %in% c('Se_Mi', 'Mod_Mi')){
    Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene, reverse = T)  
  } else {
    Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene, reverse = F)
  }
  
  DEGType <- getDEGType(Res)
  sink(paste0(pdir, '/DEGType_table.txt'))
  table(DEGType)
  sink()
  
  clu <- clusterGene(Res, gene = names(DEGType)[!DEGType %in% c('nonDEG', 'meanSig')], type = 'variable', k=5)
  # clu2 <- rep(6, sum(DEGType == 'meanSig'))
  # if (length(clu2) > 0) {
  #   names(clu2) <- names(DEGType)[DEGType %in% c('meanSig')]
  #   clu <- c(clu, clu2)
  # }
  design = Res$design
  cellanno = Res$cellanno
  meandiff <- sapply(c(0,1), function(i){
    as <- rownames(design[design[,2]==i, ])
    rowMeans(Res$expr.ori[names(DEGType)[DEGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1]])
  })
  
  large0 <- rownames(meandiff)[meandiff[,1] >= meandiff[,2]]
  large1 <- rownames(meandiff)[meandiff[,1] < meandiff[,2]]
  
  clu2 <- rep(6, length(large0))
  names(clu2) <- large0
  clu3 <- rep(7, length(large1))
  names(clu3) <- large1
  clu = c(clu, clu2, clu3)
  Res$cluster <- clu
  saveRDS(Res, paste0(rdir, paste0('numeric_res_with_clu.rds')))
  ## --------------
  ## save diff gene
  ## --------------
  allg <- diffgene
  res <- data.frame(gene = allg, statistics[allg, ], cluster = Res$cluster[allg], stringsAsFactors = F)
  res <- res[order(res[, grep('^fdr.*overall$', colnames(res))]), ]
  res <- cbind(res, DEGType = DEGType[rownames(res)])
  write.csv(res, paste0(pdir, 'differential_genes.csv'))
  
  ## ----------------
  ## plotClusterDiff
  ## ----------------
  pdf(paste0(pdir, '/cluster_diff.pdf'), width = 3, height = 2)
  plotClusterDiff(testobj=Res, gene = allg)
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
  saveRDS(goRes, paste0(pdir, '/goRes.rds'))
  
  nn <- sapply(1:length(goRes), function(i){
    tmp <- goRes[[i]]
    tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
    write.csv(tmp, paste0(pdir, 'cluster', i, '_GO.csv'))
    print(str(tmp))
    return(0)
  })
  
  pdf(paste0(pdir, '/hm_GO_term.pdf'), width = 7.2, height = 3.5)
  print(plotGOEnrich(goRes))
  dev.off()
  
  ## -----------------------
  ## plotClusterMeanAndDiff
  ## -----------------------
  pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 3.8, height = 7.5)
  print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
  dev.off()
  
  # --------------------------------------
  # compare original and fitted expression
  # --------------------------------------
  # png(paste0(pdir, 'fitHm.png'),width = 4000,height = 2200,res = 300)
  # plotFitHm(Res, type = 'variable')
  # dev.off()
  # 
  # png(paste0(pdir, 'fitHm_rownames.png'),width = 12000,height = 10000,res = 300)
  # print(plotFitHm(Res, type='variable', showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
  # dev.off()
  
  
  png(paste0(pdir, 'DiffFitHm.png'),width = 4000,height = 2200,res = 200)
  plotDiffFitHm(Res, type = 'variable', cellWidthTotal = 200, cellHeightTotal = 300)
  dev.off()
  
  # png(paste0(pdir, 'DiffFitHm_rownames.png'),width = 12000,height = 10000,res = 300)
  # print(plotDiffFitHm(Res, type='variable', showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
  # dev.off()
  
  ## ----------
  ## plot DEG 
  ## ----------
  # DEGType <- DEGType[diffgene]
  # id <- sort(sample(1:ncol(Res$populationFit[[1]]), ncol(Res$expr.ori)))
  # Res$populationFit[[1]] <- Res$populationFit[[1]][, id]
  # Res$populationFit[[2]] <- Res$populationFit[[2]][, id]
  # 
  # for (i in unique(DEGType)){  ## debug -- ok!!
  #   print(i)
  #   gene <- names(DEGType)[DEGType == i]
  #   png(paste0(pdir, 'diffgene_sampleFit_', i, '.png'), width = 4000, height = 2500, res = 200)
  #   print(plotGene(Res, gene = gene[1:min(length(gene), 25)], plot.point = T, point.size = 0.1, variable = 'type'))
  #   dev.off()
  # }
  # 
  # for (i in unique(DEGType)){
  #   print(i)
  #   gene <- names(DEGType)[DEGType == i]
  #   png(paste0(pdir, 'diffgene_groupFit_', i, '.png'), width = 2500, height = 2500, res = 200)
  #   print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)], subSampleNumber=1000))
  #   dev.off()
  # } 
  # 
  # for (i in 1:max(Res$cluster)){
  #   print(i)
  #   gene <- rownames(res)[res$cluster == i]
  #   png(paste0(pdir, 'diffgene_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
  #   print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)], subSampleNumber=1000))
  #   dev.off()
  # }
  
}


  
