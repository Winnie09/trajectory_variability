library(here)
setwd(here())
source('function/01_function.R')
ddir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
rdir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/HD_Se/'
pdir <- 'covid/Su_2020_Cell/testvar/useSaverImp/plot/HD_Se/'
dir.create(rdir, recursive = T)
dir.create(pdir, recursive = T)
Res <- readRDS(paste0(ddir, 'numeric_HD_Se_res.rds'))
Res$populationFit <- getPopulationFit(Res, gene = names(Res$fdr[Res$fdr < 0.05]), type = 'variable')
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = names(Res$fdr[Res$fdr < 0.05]))

Res$cluster <- clusterGene(Res, gene = names(Res$fdr[Res$fdr < 0.05]), type = 'variable', k=5)

allg <- names(Res$fdr[Res$fdr < 0.05])
res <- data.frame(gene = allg, pvalue = Res$pvalue[allg], fdr = Res$fdr[allg], foldchange = Res$foldchange[allg], cluster = Res$cluster[allg])
res <- res[order(res$fdr, abs(res$foldchange)), ]
write.csv(res, paste0(rdir, 'differential_genes.csv'))

goRes <- GOEnrich(testobj = Res, type = 'variable')

nn <- sapply(1:length(goRes), function(i){
  tmp <- goRes[[i]]
  tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  write.csv(tmp, paste0(rdir, 'cluster', i, '_GO.csv'))
  print(str(tmp))
  return(0)
})
  
pdf(paste0(pdir, 'cluster_mean_and_diff'), width = 5, height = 9)  
plotClusterMeanAndDiff(Res, cluster = Res$cluster)
dev.off()

for (i in 1:max(clu)){
  print(i)
  gene <- names(clu)[clu == i]
  tmp <- res[gene, ]
  png(paste0(pdir, 'diffgene_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
  plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)])
  dev.off()
}
