rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
test.method = 'EM_pm'
# for (comparison in c('HD_Se', 'HD_Mi','Mod_Se','HD_Mod', 'Se_Mi', 'Recovered_Deceased')){
# for (comparison in list.files(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method))){
for (comparison in c('Mod_Mi','Se_Mi')){ #'Recovered_Deceased', 'HD_Mi', 'HD_Se', 'Mod_Se')){ ## 'Se_Mi', 'Mod_Mi',  
  print(comparison)
  rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', test.method, '/', comparison, '/')
  pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', test.method, '/', comparison, '/')
  
  dir.create(pdir, recursive = T, showWarnings = F)
  Res <- readRDS(paste0(rdir, paste0('numeric_res.rds')))
  print(names(Res))
  
  statistics = Res$statistics
  diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
  str(diffgene)
  
  ## --------------
  ## population fit
  ## --------------
  Res$populationFit <- getPopulationFit(Res, gene = rownames(Res$expr), type = 'variable')
  
  ## -----------
  ## clustering
  ## -----------
  if (comparison %in% c('Se_Mi', 'Mod_Mi')){
    Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = rownames(Res$expr), reverse = T)  
  } else {
    Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = rownames(Res$expr), reverse = F)
  }
  
  DDGType <- getDDGType(Res)
  Res$DDGType <- DDGType
  
  ## autoclu
  # clu <- clusterGene(Res, gene = names(DDGType)[!DDGType %in% c('nonDDG')], type = 'variable', scale.difference = F, method = 'kmeans', k.auto = TRUE)
  clu <- clusterGene(Res, gene = names(Res$DDGType)[!Res$DDGType %in% c('nonDDG')], type = 'variable', scale.difference = T, method = 'kmeans', k.auto = F, k = 3) 
  sink(paste0(pdir, '/DDGType_table.txt'))
  print(table(DDGType))
  table(clu)
  sink()
  
  # 47  54 143 161 207 454 249 
  Res$cluster = clu
  saveRDS(Res, paste0(rdir, paste0('numeric_res_with_clu.rds')))
  
  ## --------------
  ## save diff gene
  ## --------------
  gd <- apply(Res$covariateGroupDiff,1, max) - apply(Res$covariateGroupDiff, 1, min)
  allg <- diffgene
  res <- data.frame(gene = allg, statistics[allg, ],cluster = Res$cluster[allg], 
                    effect_size = gd[allg], stringsAsFactors = F)
  res <- res[order(res[,2], -res[,4]), ]
  res <- cbind(res, DDGType = DDGType[rownames(res)])
  write.csv(res, paste0(pdir, 'differential_genes.csv'))
  
  ## -----------------------
  ## plotClusterMeanAndDiff
  ## -----------------------
  pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 3.8, height = 7.5)
  print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
  dev.off()
  
  ## -----------
  ## GO analysis
  ## -----------
  goRes <- GOEnrich(testobj = Res, type = 'variable')  
  saveRDS(goRes, paste0(pdir, '/goRes.rds'))
  
  nn <- sapply(names(goRes), function(i){
    tmp <- goRes[[i]]
        tmp <- tmp[tmp[, 'FDR'] < 0.25, ]
    write.csv(tmp, paste0(pdir, 'cluster', i, '_GO_FDR0.25.csv'))
    tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
    write.csv(tmp, paste0(pdir, 'cluster', i, '_GO_FDR0.05.csv'))
    print(str(tmp))
    return(0)
  })
  
  pdf(paste0(pdir, '/hm_GO_term5.pdf'), width = 6.8, height = 3.5)
  print(plotGOEnrich(goRes))
  dev.off()
  
  pdf(paste0(pdir, '/hm_GO_term10.pdf'), width = 6.8, height = 7)
  print(plotGOEnrich(goRes, n = 10))
  dev.off()
  
  pdf(paste0(pdir, '/hm_GO_term10_fdr0.25.pdf'), width = 6.8, height = 7)
  print(plotGOEnrich(goRes, n = 10, fdr.cutoff = 0.25))
  dev.off()
  
  pdf(paste0(pdir, '/hm_GO_term10_sortbyFC_fdr0.25.pdf'), width = 9, height = 7)
  print(plotGOEnrich(goRes, n = 10, sortByFDR = F,fdr.cutoff = 0.25))
  dev.off()
  
  # --------------------------------------
  # compare original and fitted expression
  # --------------------------------------
  # png(paste0(pdir, 'DiffFitHm3.png'),width = 5000,height = 3000,res = 300)
  # plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
  # dev.off()
  # 
  # png(paste0(pdir, 'DiffFitHm3_sub.png'),width = 5000,height = 3000,res = 300)
  # plotDiffFitHm3(Res,  cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = T)
  # dev.off()

  # png(paste0(pdir, 'DiffFitHm3_200.png'),width = 5000,height = 3000,res = 200)
  # plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
  # dev.off()
  # 
  # png(paste0(pdir, 'DiffFitHm3_100.png'),width = 5000,height = 3000,res = 100)
  # plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
  # dev.off()
  
  png(paste0(pdir, 'DiffFitHm5.png'),width = 10000,height = 6000,res = 300)
  plotDiffFitHm5(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
  dev.off()
  
  pdf(paste0(pdir, 'DiffFitHm5_sub1e3.pdf'),width = 30,height = 8)
  plotDiffFitHm5(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = T)
  dev.off()
  
  pdf(paste0(pdir, 'DiffFitHm5.pdf'),width = 30,height = 8)
  plotDiffFitHm5(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = F)
  dev.off()
  
  # # 
  # ########
  # png(paste0(pdir, 'DiffFitHm3_clu_changepoint.png'),width = 5000,height = 3000,res = 100)
  # plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
  # dev.off()
  # 
  # png(paste0(pdir, 'DiffFitHm3_clu_magnitute_cor.png'),width = 5000,height = 3000,res = 100)
  # plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
  # dev.off()
  # 
  # png(paste0(pdir, 'DiffFitHm3_clu_type_magnitute_changepoint.png'),width = 5000,height = 3000,res = 100)
  # plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
  # dev.off()
  # 
  
  
  groupdiff <- Res$covariateGroupDiff
  mat1 <- groupdiff[names(clu)[clu == max(clu)-1], ]
  mat2 <- groupdiff[names(clu)[clu == max(clu)], ]
  dn1 <- dimnames(mat1)
  dn2 <- dimnames(mat2)
  v1 <- as.vector(mat1)
  v2 <- as.vector(mat2)
  mat1 <- matrix(v1/max(v1), nrow(mat1))
  mat2 <- matrix(v2/max(v2), nrow(mat2))
  dimnames(mat1) <- dn1
  dimnames(mat2) <- dn2
  
  d1 <- apply(mat1, 1, max) - apply(mat1, 1, min)
  d2 <- apply(mat2, 1, max) - apply(mat2, 1, min)
  mat1 <- mat1[order(d1), ]  
  mat2 <- mat2[order(d2), ]  
  Res$covariateGroupDiff <- rbind(groupdiff[names(clu)[clu %in% seq(1, max(clu)-2)], ], mat1, mat2)
   
  # png(paste0(pdir, 'DiffFitHm3_clu_type_changepoint_cor.png'),width = 5000,height = 3000,res = 100)
  # plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = F)
  # dev.off()
  
  # png(paste0(pdir, 'DiffFitHm_rownames3.png'),width = 12000,height = 10000,res = 300)
  # print(plotDiffFitHm3(Res, type='variable', showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
  # dev.off()

  # ## ----------
  # ## plot DDG
  # ## ----------
  DDGType <- DDGType[diffgene]
  # id <- sort(sample(1:ncol(Res$populationFit[[1]]), ncol(Res$expr)))
  # Res$populationFit[[1]] <- Res$populationFit[[1]][, id]
  # Res$populationFit[[2]] <- Res$populationFit[[2]][, id]
  # colnames(Res$population[[1]]) <- colnames(Res$population[[2]]) <- id

  # for (i in unique(DDGType)){  ## debug -- ok!!
  #   print(i)
  #   gene <- names(DDGType)[DDGType == i]
  #   png(paste0(pdir, 'diffgene_sampleFit_', i, '.png'), width = 4000, height = 2500, res = 200)
  #   print(plotGene(Res, gene = gene[1:min(length(gene), 25)], plot.point = T, point.size = 0.1, variable = 'type'))
  #   dev.off()
  # }
  # 
  # for (i in unique(DDGType)){
  #   print(i)
  #   gene <- names(DDGType)[DDGType == i]
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





