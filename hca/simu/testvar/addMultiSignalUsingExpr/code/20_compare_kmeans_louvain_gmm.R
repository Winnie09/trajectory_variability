library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))
library(ComplexHeatmap)
library(mclust)
for (dataType  in seq(2,4)){
  print(dataType)
  ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr'
  rdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/cluster/',dataType)
  pdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/plot/EM_pm/', dataType)
  dir.create(pdir, recursive = T, showWarnings = F)
  dir.create(rdir, recursive = T, showWarnings = F)
  Res <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/EM_pm/',dataType,'.rds'))
  names(Res)
  stat <- Res$statistics
  head(stat)
  stat <- stat[order(stat[,1], -stat[,3]), ]
  ## identify XDE genes with FDR.overall < 0.05 cutoff
  diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
  str(diffgene)
  Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')
  Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene)
  colnames(Res$populationFit[[1]]) <- colnames(Res$populationFit[[2]]) <- colnames(Res$expr)
  clu_louvain <- clusterGene(Res, gene = diffgene, type = 'variable', k=7, method = 'louvain')
  saveRDS(clu_louvain, paste0(rdir, '/cluster_louvain.rds'))
  if (dataType == 4){
    clu_kmeans <- clusterGene(Res, gene = diffgene, type = 'variable', k=6, method = 'kmeans')
    clu_gmm <- clusterGene(Res, gene = diffgene, type = 'variable', k=6, method = 'GMM')
  } else {
    clu_kmeans <- clusterGene(Res, gene = diffgene, type = 'variable', k=7, method = 'kmeans')
    clu_gmm <- clusterGene(Res, gene = diffgene, type = 'variable', k=7, method = 'GMM')
  }
    
  saveRDS(clu_kmeans, paste0(rdir, '/cluster_kmeans.rds'))
  # clu_kmeans = readRDS(paste0(rdir, '/cluster_kmeans.rds'))
  saveRDS(clu_gmm, paste0(rdir, '/cluster_gmm.rds'))
  # clu_gmm = readRDS(paste0(rdir, '/cluster_gmm.rds'))
  
  ## ----------------
  ## plot confusion matrix
  ## ----------------
  m <- caret::confusionMatrix(data = as.factor(paste0('cluster',clu_louvain)), reference = as.factor(paste0('cluster', clu_kmeans[names(clu_louvain)])),
                              dnn = c('louvain','kmeans'))  ## dnn = c("Prediction", "Reference")
  tb <- t(m$table)
  tb <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm_louvain_kmeans.pdf'), width = 3, height = 2.57)
  print(pheatmap(tb, cluster_rows = F, cluster_cols = F))
  dev.off()

  m <- caret::confusionMatrix(data = as.factor(paste0('cluster',clu_gmm)), reference = as.factor(paste0('cluster', clu_kmeans[names(clu_gmm)])),
                              dnn = c('GMM','kmeans'))  ## dnn = c("Prediction", "Reference")
  tb <- t(m$table)
  tb <- tb/rowSums(tb)
  pdf(paste0(pdir, '/confusionMatrix_hm_gmm_kmeans.pdf'), width = 3, height = 2.57)
  print(pheatmap(tb, cluster_rows = F, cluster_cols = F))
  dev.off()

  ## -----------------------------------
  ## plot fitted gene expression heatmap
  ## -----------------------------------
  testobj = Res
  showRowName = FALSE
  cellWidthTotal = 250
  cellHeightTotal = 400
           showCluster = T
           colann = NULL
           rowann = NULL
           annotation_colors = NULL
           subsampleCell = TRUE
           numSubsampleCell = 1e3
           sep = NA
           break.0 = TRUE

  library(gridExtra)
  Res$cluster <- clu_kmeans
  png(
    paste0(pdir, '/hm_kmeans.png'),
    bg = 'white',
    width = 4500,
    height = 3500,
    res = 300
  )
  plotXDEHm(
    Res,
    subsampleCell  = F,
    showCluster = TRUE,
    cellWidthTotal = 130,
    cellHeightTotal = 200
  )
  dev.off()

  Res$cluster <- clu_louvain
  png(
    paste0(pdir, '/hm_louvain.png'),
    bg = 'white',
    width = 4500,
    height = 3500,
    res = 300
  )
  plotXDEHm(
    Res,
    subsampleCell  = F,
    showCluster = TRUE,
    cellWidthTotal = 130,
    cellHeightTotal = 200
  )
  dev.off()

  library(RColorBrewer)
  Res$cluster <- clu_gmm
  png(
    paste0(pdir, '/hm_gmm.png'),
    bg = 'white',
    width = 4500,
    height = 3500,
    res = 300
  )
  library(circlize)
  plotXDEHm(
    Res,
    subsampleCell  = F,
    showCluster = TRUE,
    cellWidthTotal = 130,
    cellHeightTotal = 200
  )
  dev.off()
}


