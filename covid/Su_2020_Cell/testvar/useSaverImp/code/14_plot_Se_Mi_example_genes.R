rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
ddir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'

comparison = 'Se_Mi'
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

res <- res[names(Res$cluster), ]
res <- cbind(res, cluster = Res$cluster)
  
## group fit    
png(paste0(pdir, 'diffgene_groupFit_example_cluster1.png'), width = 1000, height = 800, res = 200)
print(plotGenePopulation(testobj = Res, type = 'variable', gene = c('IFI35', 'IRF7', 'MX2')))
dev.off()

png(paste0(pdir, 'diffgene_point_example_cluster1.png'), width = 1400, height = 500, res = 200)
plotGene(testptObj = Res, variable = 'type', gene = c('IFI35', 'IRF7', 'MX2'), free.scale = TRUE, facet.sample = FALSE, plot.point = T, line.alpha = 1, line.size = 0.5, point.alpha=0.2, point.size=0.2, original.expr = TRUE)
dev.off()

png(paste0(pdir, 'diffgene_groupDiff_example_cluster1.png'),width = 2500,height = 2500, res = 200)
print(plotClusterDiff(testobj = Res, gene = c('IFI35', 'IRF7', 'MX2'), each = TRUE, sep = ':.*'))
dev.off()

## cluster 2
png(paste0(pdir, 'diffgene_groupFit_example_cluster2.png'), width = 1000, height = 800, res = 200)
print(plotGenePopulation(testobj = Res, type = 'variable', gene = c('FLI1', 'RELB', 'IL15RA')))
dev.off()

png(paste0(pdir, 'diffgene_groupDiff_example_cluster2.png'),width = 2500,height = 2500, res = 200)
print(plotClusterDiff(testobj = Res, gene = c('FLI1', 'RELB', 'IL15RA'), each = TRUE, sep = ':.*'))
dev.off()

png(paste0(pdir, 'diffgene_point_example_cluster2.png'), width = 1400, height = 500, res = 200)
plotGene(testptObj = Res, variable = 'type', gene = c('FLI1', 'RELB', 'IL15RA'), free.scale = TRUE, facet.sample = FALSE, plot.point = T, line.alpha = 1, line.size = 0.5, point.alpha=0.2, point.size=0.2, original.expr = TRUE)
dev.off()

## cluster 3
png(paste0(pdir, 'diffgene_groupFit_example_cluster3.png'), width = 1000, height = 800, res = 200)
print(plotGenePopulation(testobj = Res, type = 'variable', gene = c('PSPH', 'YPEL3', 'CD37')))
dev.off()  

png(paste0(pdir, 'diffgene_groupDiff_example_cluster3.png'),width = 2500,height = 2500, res = 200)
print(plotClusterDiff(testobj = Res, gene = c('PSPH', 'YPEL3', 'CD37'), each = TRUE, sep = ':.*'))
dev.off()

png(paste0(pdir, 'diffgene_point_example_cluster3.png'), width = 1400, height = 500, res = 200)
plotGene(testptObj = Res, variable = 'type', gene = c('PSPH', 'YPEL3', 'CD37'), free.scale = TRUE, facet.sample = FALSE, plot.point = T, line.alpha = 1, line.size = 0.5, point.alpha=0.2, point.size=0.2, original.expr = TRUE)
dev.off()

