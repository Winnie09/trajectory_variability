rm(list=ls())
library(here)
setwd(here())
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/res/twogroup/'
i = as.numeric(commandArgs(trailingOnly = T)[[1]])
comp = sub('.rds', '', list.files(ddir, pattern = '.rds'))[i]
## "Frontal_Occipital": frontal 1, occipital 0 (higher expression)
pdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/brain/atlasGBM/GBMonly/trajectory/plot/twogroup/', comp, '/')
dir.create(pdir, showWarnings = F, recursive = T)

Res <- readRDS(paste0(ddir, comp, '.rds'))
print(names(Res))
str(Res$llr.overall)
stat = Res$statistics
diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
str(diffgene)
stat <- stat[order(stat[,1], -stat[,3]), ]
## if the above test works well, then refine the following codes
## --------------
## population fit
## --------------
Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')

## -----------
## clustering
## -----------
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene, reverse = T)

DDGType <- getDDGType(Res)
print(table(DDGType))
Res$DDGType <- DDGType
sink(paste0(pdir, 'DDGType_table.txt'))
print(table(DDGType))
sink()

clu <- clusterGene(Res, gene = names(Res$DDGType)[!Res$DDGType %in% c('nonDDG')], type = 'variable', k.auto=F, scale.difference = T, k = 2)
table(clu)
Res$cluster = clu
# saveRDS(Res, paste0(pdir, paste0('numeric_res_with_clu.rds')))
## --------------
## save diff gene
## --------------
allg <- diffgene
res <- data.frame(gene = allg, stat[allg, ], cluster = Res$cluster[allg], stringsAsFactors = F)
res <- res[order(res[,2], -res[,4]), ]
res <- cbind(res, DDGType = DDGType[rownames(res)])
write.csv(res, paste0(pdir, 'differential_genes.csv'))

# ## -----------------------
# ## plotClusterMeanAndDiff
# ## -----------------------
# pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 3.8, height = 7.5)
# print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
# dev.off()

## -----------
## GO analysis
## -----------
goRes <- GOEnrich(testobj = Res, type = 'variable')  
# saveRDS(goRes, paste0(pdir, '/goRes.rds'))

nn <- sapply(names(goRes), function(i){
  tmp <- goRes[[i]]
  tmp <- tmp[tmp[, 'FDR'] < 0.05, , drop = FALSE]
  write.csv(tmp, paste0(pdir, 'cluster', i, '_GO.csv'))
  print(str(tmp))
  return(nrow(tmp))
})

if (sum(nn) > 0){
  pdf(paste0(pdir, 'GO.pdf'), width = 8, height = 5)
  plotGOEnrich(goRes)
  dev.off()
}

png(paste0(pdir, 'DiffFitHm5_rownames.png'),width = 12000,height = 3500,res = 300)
plotDiffFitHm5(Res, showRowName = T, cellWidthTotal = 230, cellHeightTotal = length(Res$cluster) * 10, sep = ':.*' , subsampleCell = FALSE, break.0 = FALSE)
dev.off()  

pdf(paste0(pdir, 'DiffFitHm5_rownames.pdf'),width = 50,height = 12)
plotDiffFitHm5(Res, showRowName = T, cellWidthTotal = 230, cellHeightTotal = length(Res$cluster) * 10, sep = ':.*' , subsampleCell = FALSE, break.0 = FALSE)
dev.off()


