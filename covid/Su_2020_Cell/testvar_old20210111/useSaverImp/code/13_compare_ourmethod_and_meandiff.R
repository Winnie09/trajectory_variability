library(here)
setwd(here())
source('function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
for (comparison in c('HD_Se', 'HD_Mi','Mod_Se','HD_Mod', 'Se_Mi', 'Recovered_Deceased')){
  pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', comparison, '/')
  rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', comparison, '/')
  res.lm <- readRDS(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/numeric_', comparison, '_meandiff.rds'))
  res.lm <- res.lm[res.lm[,5] < 0.05, ]
  Res <- readRDS(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/numeric_', comparison, '_res.rds'))
  res <- data.frame(fdr = Res$fdr, foldchange = Res$foldchange, stringsAsFactors = FALSE)
  rownames(res) <- names(Res$fdr)
  res <- res[res[,1] < 0.05, ]
  res <- res[order(res[,1], -res[,2]), ]
  g.new <- rownames(res)[!rownames(res) %in% rownames(res.lm)]
  
  allg <- rownames(res)
  o1 <- seq(1,nrow(res))
  res.lm.tmp <- readRDS(paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/numeric_', comparison, '_meandiff.rds'))
  o2 <- unlist(sapply(allg, function(g) which(rownames(res.lm.tmp) == g)))
  pdf(paste0(pdir, 'diffgene_order_compared_with_limma_order.pdf'), width = 4, height = 3.8)
  plot(o1~o2, xlab = 'limma order', ylab = 'ourmethod order', pch = 20, main='differential genes')
  abline(v = nrow(res.lm), col='red')
  dev.off()
  
  tb <- read.csv(paste0(rdir, 'differential_genes.csv'), row.names = 1, as.is = TRUE)
  table(tb[tb[,1] %in% g.new, 5])
  tb.new <- cbind(tb, limma = ifelse(tb[,1] %in% rownames(res.lm), 'meandiff', 'trenddiff'))
  write.csv(tb.new, paste0(rdir, 'differential_genes_marked_with_limma.csv'))
  
  ourClu1 <- rownames(tb[tb[,5] == 1, ])
  ourClu2 <- rownames(tb[tb[,5] == 2, ])
  ourClu3 <- rownames(tb[tb[,5] == 3, ])
  meandiff <- rownames(res.lm)
  
  library(VennDiagram)
  venn.diagram(
    x = list(ourClu1, ourClu2, ourClu3, meandiff),
    category.names = c("ourClu1" , "ourClu2 " , "ourClu3", 'meandiff'),
    filename = paste0(pdir, 'venn_diagramm.png'),
    output=TRUE
  )
  
  ######## pipeline
  ## --------------
  ## population fit
  ## --------------
  Res$populationFit <- getPopulationFit(Res, gene = g.new, type = 'variable')
  
  ## -----------
  ## clustering
  ## -----------
  Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = g.new)
  clu <- tb[g.new, 5]
  names(clu) <- g.new
  Res$cluster <- clu
  
  ## --------------
  ## save diff gene
  ## --------------
  write.csv(tb[g.new, ], paste0(rdir, 'differential_genes_not_detected_by_limma.csv'))
  
  ## ----------------
  ## plotClusterDiff
  ## ----------------
  pdf(paste0(pdir, 'diffgeneNotDetectedByLimma_cluster_diff.pdf'), width = 3, height = 2)
  plotClusterDiff(testobj=Res, gene = g.new)
  dev.off()
  
  ## ---------------
  ## plotClusterMean
  ## ----------------
  pdf(paste0(pdir, 'diffgeneNotDetectedByLimma_cluster_mean.pdf'), width = 5, height = 3.5)
  plotClusterMean(testobj=Res, cluster = Res$cluster)
  dev.off()
  
  ## -----------
  ## GO analysis
  ## -----------
  goRes <- GOEnrich(testobj = Res, type = 'variable', use.clusters = FALSE)
  tmp <- goRes[[1]]
  tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  write.csv(tmp, paste0(rdir, 'diffgeneNotDetectedByLimma_GO.csv'))
  
  
  goRes <- GOEnrich(testobj = Res, type = 'variable', use.clusters = TRUE)
  for (i in 1:length(goRes)){
    tmp <- goRes[[i]]
    # tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
    tmp <- tmp[1:50, ]
    write.csv(tmp, paste0(rdir, 'diffgeneNotDetectedByLimma_GO_cluster', i, '.csv'))
  }
  
  ## -----------------------
  ## plotClusterMeanAndDiff
  ## -----------------------
  pdf(paste0(pdir, 'diffgeneNotDetectedByLimma_cluster_mean_and_diff.pdf'), width = 4, height = 5.5)
  print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
  dev.off()
  
  for (i in 1:max(Res$cluster)){
    print(i)
    gene <- g.new
    png(paste0(pdir, 'diffgeneNotDetectedByLimma_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
    print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
    dev.off()
  }
  
  for (i in 1:max(Res$cluster)){
    print(i)
    gene <- g.new
    png(paste0(pdir, 'diffgeneNotDetectedByLimma__groupDiff_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
    print(plotClusterDiff(testobj = Res, gene = gene[1:min(length(gene), 100)], each = TRUE, sep = ':.*'))
    dev.off()
  }
}
  

