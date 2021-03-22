rm(list = ls())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
library(here)
setwd(here())
source('function/01_function.R')
Res <- readRDS('covid/Su_2020_Cell/testtime/result/EM_pm/testtime_res.rds')
pdir <- rdir <- 'covid/Su_2020_Cell/testtime/plot/EM_pm/'
dir.create(pdir, recursive = T, showWarnings = F)
### downstream 
stat = Res$statistics 
str(stat)
str(Res$expr)

pdir <- 'covid/Su_2020_Cell/testtime/plot/EM_pm/'
# ---------------------------- #
# downstream analysis pipeline #  
# ---------------------------- #
diffgene = rownames(Res$statistics)[Res$statistics[,grep('^fdr.*overall$', colnames(Res$statistics))] < 0.05]
str(diffgene)

sum(stat[,1]< 0.05)
## --------------
## population fit
## --------------
Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'time')

## -----------
## clustering
## -----------
Res$cluster <- clusterGene(Res, gene = diffgene, type = 'time', k=5)
saveRDS(Res, 'covid/Su_2020_Cell/testtime/result/EM_pm/testtime_res_with_cluster.rds')
## --------------
## save diff gene
## --------------
allg <- rownames(Res$statistics[Res$statistics[,1]<0.05,,drop=FALSE])
res <- Res$statistics[allg, ]
res <- res[order(res[,1], -res[,3]), ]
res <- cbind(res, cluster = Res$cluster[rownames(res)])
write.csv(res, paste0(pdir, 'testtime_differential_genes.csv'))

## ---------------
## plotClusterMean
## ----------------
pdf(paste0(pdir, 'cluster_mean.pdf'), width = 5, height = 3.5)
plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'time')
dev.off()

## -----------
## GO analysis
## -----------
goRes <- GOEnrich(testobj = Res, type = 'time', sep = ':.*')
saveRDS(goRes, paste0(pdir, 'goRes.rds'))

nn <- sapply(1:length(goRes), function(i){
  tmp <- goRes[[i]]
  # tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  write.csv(tmp, paste0(pdir, 'cluster', i, '_GO.csv'))
  saveRDS(tmp, paste0(pdir, 'cluster', i, '_GO.rds'))
  print(str(tmp))
  return(0)
})

pdf(paste0(pdir, 'hm_GO_term.pdf'), width = 7.2, height = 3.5)
print(plotGOEnrich(goRes))
dev.off()

# ------------------------------------------------------
# compare original and fitted expression: not tested yet
# ------------------------------------------------------
colnames(Res$populationFit) <- colnames(Res$expr)
png(paste0(pdir, 'fitHm.png'),width = 3000,height = 2500,res = 200)
plotFitHm(Res, cellHeightTotal = 250, cellWidthTotal=200, subsampleCell=FALSE )
dev.off()

png(paste0(pdir, 'fitHm_rownames.png'),width = 12000,height = 10000,res = 200)
print(plotFitHm(Res, showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
dev.off()









