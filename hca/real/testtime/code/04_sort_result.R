library(here)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(grid)
setwd(here())
path = 'lymph'
rdir <- ddir <- paste0('hca/real/testtime/result/', path)
pdir  <- paste0('hca/real/testtime/plot/', path)
dir.create(pdir, recursive = T)
source('function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
Res <- readRDS(paste0(ddir, '/testtime_res.rds'))

## get the differential genes
res = Res$statistics
res <- res[order(res[,1], -res[,2]), ]
diffgene <- rownames(res)[res[,1] < 0.05]

## --------------
## population fit
## --------------
# fit <- t(sapply(rownames(res)[res$fdr < 0.05], function(g){
#   tmp <- get_population_fit(Res, variable = NA, gene = g)
# }))

Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'time')
saveRDS(Res$populationFit, paste0(rdir, '/population_fit.rds')) 

## -----------
## clustering
## -----------
Res$cluster <- clusterGene(Res, gene = diffgene, type = 'time', k=5)
saveRDS(Res$cluster, paste0(rdir, '/cluster.rds'))
  
## --------------
## save diff gene
## --------------
allg <- rownames(Res$statistics[Res$statistics[,1]<0.05,,drop=FALSE])
res <- Res$statistics[allg, ]
res <- res[order(res$fdr, -res$fc), ]
write.csv(res, paste0(rdir, '/testtime_differential_genes.csv'))

## ---------------
## plotClusterMean
## ----------------
pdf(paste0(pdir, '/cluster_mean.pdf'), width = 5, height = 3.5)
plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'time')
dev.off()

## -----------
## GO analysis
## -----------
goRes <- GOEnrich(testobj = Res, type = 'time', version = 1, sep = ':.*')
saveRDS(goRes, paste0(rdir, '/goRes.rds'))

nn <- sapply(1:length(goRes), function(i){
  tmp <- goRes[[i]]
  # tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  write.csv(tmp, paste0(rdir, '/cluster', i, '_GO.csv'))
  saveRDS(tmp, paste0(rdir, '/cluster', i, '_GO.rds'))
  print(str(tmp))
  return(0)
})

pdf(paste0(pdir, '/hm_GO_term.pdf'), width = 7.2, height = 3.5)
print(plotGOEnrich(goRes))
dev.off()

# ------------------------------------------------------
# compare original and fitted expression: not tested yet
# ------------------------------------------------------
png(paste0(pdir, '/fitHm.png'),width = 3000,height = 2500,res = 300)
print(plotFitHm(Res))
dev.off()

png(paste0(pdir, '/fitHm_rownames.png'),width = 12000,height = 10000,res = 300)
print(plotFitHm(Res, showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
dev.off()



