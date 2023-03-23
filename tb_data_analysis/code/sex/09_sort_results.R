rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')

rdir <- paste0('tb/res/sex/pc2/')
pdir <- paste0('tb/plot/sex/')
dir.create(pdir, recursive = T, showWarnings = F)
Res <- readRDS('tb/res/sex/pc2/lamian_pm_allcores.rds')

statistics = Res$statistics
diffgene <- rownames(statistics[statistics[, grep('^fdr.*overall$', colnames(statistics))] < 0.05,])
str(diffgene)

## --------------
## population fit
## --------------
Res$expr <- readRDS('tb/data/sex/expr.rds')
Res$testvar = ncol(Res$design)
# Res$design = Res$design[, c(1, 39, 2:38)]
Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')
saveRDS(Res$populationFit, 'tb/res/sex/pc2/lamian_pm_populationFit.rds')

## -----------
## clustering
## -----------
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene, reverse = T)  

XDEType <- getXDEType(Res)
Res$XDEType <- XDEType

## autoclu
# clu <- clusterGene(Res, gene = names(XDEType)[!XDEType %in% c('nonXDE')], type = 'variable', scale.difference = F, method = 'kmeans', k.auto = TRUE)
clu <- clusterGene(Res, gene = names(Res$XDEType)[!Res$XDEType %in% c('nonXDE')], type = 'variable', scale.difference = T, method = 'kmeans', k.auto = F, k = 3) 
sink(paste0(pdir, '/XDEType_table.txt'))
print(table(XDEType))
table(clu)
sink()

# 47  54 143 161 207 454 249 
Res$cluster = clu

## --------------
## save diff gene
## --------------
gd <- apply(Res$covariateGroupDiff,1, max) - apply(Res$covariateGroupDiff, 1, min)
allg <- diffgene
res <- data.frame(gene = allg, statistics[allg, ],cluster = Res$cluster[allg], 
                  effect_size = gd[allg], stringsAsFactors = F)
res <- res[order(res[,2], -res[,4]), ]
res <- cbind(res, XDEType = XDEType[rownames(res)])
write.csv(res, paste0(pdir, 'differential_genes.csv'))

## -----------------------
## plotClusterMeanAndDiff
## -----------------------
# pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 3.8, height = 7.5)
# print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
# dev.off()

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
png(paste0(pdir, 'DiffFitHm5.png'),width = 2500,height = 2200,res = 300)
plotDiffFitHm5(Res, subsampleCell = TRUE)
dev.off()

Res$expr <- NULL
saveRDS(Res, paste0(rdir, paste0('numeric_res_with_clu.rds')))

id = seq(from = 1, to = length(Res$pseudotime), length.out = ncol(Res$populationFit[[1]]))
colnames(Res$populationFit[[1]]) <- colnames(Res$populationFit[[2]]) <- colnames(Res$expr)[id]

pdf(paste0(pdir, 'cluster_mean.pdf'), width = 5, height = 3.5)
plotClusterMean(testobj = Res, type = 'Variable', facet = TRUE)
dev.off()

######################################################################
## plot cluster mean
## split up cluster to cluster a and b (e.g. cluster 1 to 1a, and 1b)
######################################################################
abtype = unlist(sapply(Res$XDEType[names(Res$cluster)], function(i){
  if (i == 'trendSig') {
    'a'
  } else if (i == 'bothSig'){
    'b'
  } else {
    ''
  }
}))
clutype = paste0(Res$cluster, abtype)
names(clutype) <- names(Res$cluster)
clutype = clutype[!clutype %in% seq(1,6)]

pdf(paste0(pdir, 'cluster_mean_abtype.pdf'), width = 3.8, height = 6.8)
plotClusterMean(testobj = Res, cluster = clutype, type = 'Variable', facet = TRUE, facet_scales = 'free', facet_nrow = 7)
dev.off()

pdf(paste0(pdir, 'cluster_group_difference.pdf'), width = 2.2, height = 7.5)
plotClusterDiff(testobj = Res, 
                            gene = names(Res$cluster),
                            cluster = Res$cluster,
                            each = T,
                            facet_scales = 'fixed',
                            facet_variable = 'cluster',
                            facet_nrow = 8,
                            sep = '',
                            reverse = F)
dev.off()

