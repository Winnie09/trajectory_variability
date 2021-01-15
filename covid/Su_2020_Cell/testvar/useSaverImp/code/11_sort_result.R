library(here)
setwd(here())
source('function/01_function.R')
ddir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
# for (comparison in c('HD_Se', 'HD_Mi','Mod_Se','HD_Mod', 'Se_Mi', 'Recovered_Deceased')){
for (comparison in c('Se_Mi', 'Recovered_Deceased')){
  print(comparison)
  rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', comparison, '/')
  pdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/', comparison, '/')
  dir.create(rdir, recursive = T)
  dir.create(pdir, recursive = T)
  Res <- readRDS(paste0(ddir, paste0('numeric_', comparison, '_res.rds')))
  
  ## --------------
  ## population fit
  ## --------------
  Res$populationFit <- getPopulationFit(Res, gene = names(Res$fdr[Res$fdr < 0.05]), type = 'variable')
  
  ## -----------
  ## clustering
  ## -----------
  Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = names(Res$fdr[Res$fdr < 0.05]))
  Res$cluster <- clusterGene(Res, gene = names(Res$fdr[Res$fdr < 0.05]), type = 'variable', k=3)
  
  ## --------------
  ## save diff gene
  ## --------------
  allg <- names(Res$fdr[Res$fdr < 0.05])
  res <- data.frame(gene = allg, pvalue = Res$pvalue[allg], fdr = Res$fdr[allg], foldchange = Res$foldchange[allg], cluster = Res$cluster[allg])
  res <- res[order(res$fdr, abs(res$foldchange)), ]
  write.csv(res, paste0(rdir, 'differential_genes.csv'))
  
  ## ----------------
  ## plotClusterDiff
  ## ----------------
  pdf(paste0(pdir, '/cluster_diff.pdf'), width = 3, height = 2)
  plotClusterDiff(testobj=Res, gene = names(Res$fdr[Res$fdr < 0.05]))
  dev.off()
  
  ## ---------------
  ## plotClusterMean
  ## ----------------
  pdf(paste0(pdir, '/cluster_mean.pdf'), width = 5, height = 3.5)
  plotClusterMean(testobj=Res, cluster = Res$cluster)
  dev.off()
  
  ## -----------
  ## GO analysis
  ## -----------
  goRes <- GOEnrich(testobj = Res, type = 'variable')
  nn <- sapply(1:length(goRes), function(i){
    tmp <- goRes[[i]]
    tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
    write.csv(tmp, paste0(rdir, 'cluster', i, '_GO.csv'))
    print(str(tmp))
    return(0)
  })
  
  ## -----------------------
  ## plotClusterMeanAndDiff
  ## -----------------------
  pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 4, height = 5.5)
  print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
  dev.off()
  
  for (i in 1:max(Res$cluster)){
    print(i)
    gene <- rownames(res[res$cluster == i, ])
    png(paste0(pdir, 'diffgene_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
    print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)]))
    dev.off()
  }

  for (i in 1:max(Res$cluster)) {
    print(i)
    gene <- rownames(res[res$cluster == i,])
    png(paste0(pdir, 'diffgene_groupDiff_cluster', i, '.png'),width = 2500,height = 2500, res = 200)
    print(plotClusterDiff(testobj = Res, gene = gene[1:min(length(gene), 100)], each = TRUE, sep = ':.*'))
    dev.off()
  }

  # --------------------------------------
  # compare original and fitted expression
  # --------------------------------------
  png(paste0(pdir, 'fitHm.png'),width = 3000,height = 2500,res = 300)
  print(plotFitHm(Res))
  dev.off()

  png(paste0(pdir, 'fitHm_rownames.png'),width = 11080,height = length(Res$cluster)*25,res = 200)
  print(plotFitHm(Res, showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 8))
  dev.off()
}
  
