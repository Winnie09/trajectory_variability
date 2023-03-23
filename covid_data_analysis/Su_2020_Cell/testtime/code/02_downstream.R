rm(list = ls())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
library(here)
setwd(here())
source('function/01_function.R')
pdir <- 'covid/Su_2020_Cell/testtime/plot/EM_pm/'
Res <- readRDS('covid/Su_2020_Cell/testtime/result/EM_pm/testtime_res.rds')
pdir <- rdir <- 'covid/Su_2020_Cell/testtime/plot/EM_pm/'
dir.create(pdir, recursive = T, showWarnings = F)
### downstream 
stat = Res$statistics 
str(stat)
str(Res$expr)


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

## plot example genes

tb = read.csv('/Users/wenpinhou/Dropbox/trajectory_variability/covid/Su_2020_Cell/testtime/plot/EM_pm/marked_DEG_red_chen.csv')
str(tb)
examplegene = tb[,1]
examplegene = names(sort(Res$cluster[tb[,1]]))

png(paste0(pdir, 'examplegene_sampleFit_all.png'), res = 200, width = 4000, height = 2000)
plotGene(Res, examplegene, plot.point = T, point.size = 0.05, point.alpha = 0.5, line.size = 0.5, line.alpha = 0.8, axis.text.blank = F)
dev.off()

png(paste0(pdir, 'examplegene_populationFit_all.png'), res = 100, width = 1000, height = 1000)
plotGenePopulation(Res, examplegene, facet.grid = T, axis.text.blank = T)
dev.off()


# gene <- c('TCF7', 'SELL', 'JUNB', 'NFKBIA', 'CD7', 'GZMA', 'CCL5', 'GZMB', 'GNLY')
# png(paste0(pdir, 'examplegene_sampleFit.png'), res = 300, width = 4500, height = 1300)
# plotGene(Res, gene, plot.point = T, point.size = 0.05, point.alpha = 0.5, line.size = 0.5, line.alpha = 0.8, axis.text.blank = F, free.scale = F)
# dev.off()
# 

# png(paste0(pdir, 'examplegene_populationFit.png'), res = 300, width = 1300, height = 1000)
# plotGenePopulation(Res, gene, axis.text.blank = F, free.scale = F, line.size = 2)
# dev.off()


# gene <- c('TCF7', 'SELL', 'JUNB', 'CD7', 'GZMA', 'CCL5', 'GZMB', 'GNLY')
gene <- c('TCF7', 'SELL', 'JUNB', 'CD7', 'CCL5', 'IFNG','GZMB', 'GNLY')
png(paste0(pdir, 'examplegene_sampleFit2.png'), res = 300, width = 3000, height = 1300)
plotGene(Res, gene, plot.point = T, point.size = 0.05, point.alpha = 0.5, line.size = 0.1, line.alpha = 0.8, axis.text.blank = F, free.scale = F, ncol = 2, legend.position = 'right')
dev.off()

pdf(paste0(pdir, 'examplegene_sampleFit2.pdf'), width = 12, height = 4.8)
plotGene(Res, gene, plot.point = T, point.size = 0.01, point.alpha = 0.05, line.size = 0.1, line.alpha = 0.7, axis.text.blank = F, free.scale = F, ncol = 2, legend.position = 'right')
dev.off()


png(paste0(pdir, 'examplegene_populationFit2.png'), res = 300, width = 2000, height = 800)
plotGenePopulation(Res, gene, axis.text.blank = F, free.scale = F, line.size = 2, ncol = 4)
dev.off()



