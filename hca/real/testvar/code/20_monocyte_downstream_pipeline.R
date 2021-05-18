rm(list=ls())
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
path = 'monocyte'
pdir <- paste0('hca/real/testvar/plot/EM_pm/', path, '/gender/')
dir.create(pdir, showWarnings = F, recursive = T)
## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')

Res <- readRDS(paste0('hca/real/testvar/result/EM_pm/', path, '/gender/gender_res.rds'))
Res$statistics <- Res$statistics[rownames(Res$statistics)!='RPS4Y1:ENSG00000129824',]
print(names(Res))
stat = Res$statistics
diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
str(diffgene)
stat <- stat[order(stat[,1], -stat[,3]), ]
## if the above test works well, then refine the following codes
## --------------
## population fit
## --------------
Res$populationFit <- getPopulationFit(Res, gene = rownames(Res$expr), type = 'variable')

## -----------
## clustering
## -----------
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = rownames(Res$expr), reverse = F)

DDGType <- getDDGType(Res)
Res$DDGType <- DDGType
sink(paste0(pdir, '/DDGType_table.txt'))
print(table(DDGType))
sink()

clu <- clusterGene(Res, gene = names(Res$DDGType)[!Res$DDGType %in% c('nonDDG')], type = 'variable', k.auto=F, scale.difference = T, k = 3)
table(clu)
Res$cluster = clu
saveRDS(Res, paste0(pdir, paste0('numeric_res_with_clu.rds')))
## --------------
## save diff gene
## --------------
allg <- diffgene
res <- data.frame(gene = allg, stat[allg, ], cluster = Res$cluster[allg], stringsAsFactors = F)
res <- res[order(res[,2], -res[,4]), ]
res <- cbind(res, DDGType = DDGType[rownames(res)], chrX = ifelse(sub(':.*', '', rownames(res)) %in% u1, T, F), chrY = ifelse(sub(':.*', '', rownames(res)) %in% u2, T, F))
write.csv(res, paste0(pdir, 'differential_genes.csv'))

## -----------------------
## plotClusterMeanAndDiff
## -----------------------
pdf(paste0(pdir, 'cluster_mean_and_diff.pdf'), width = 3.8, height = 7.5)
print(plotClusterMeanAndDiff(Res, cluster = Res$cluster))
dev.off()

## -----------
## GO analysis
## -----------
goRes <- GOEnrich(testobj = Res, type = 'variable')  
saveRDS(goRes, paste0(pdir, '/goRes.rds'))

nn <- sapply(names(goRes), function(i){
  tmp <- goRes[[i]]
  tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  write.csv(tmp, paste0(pdir, 'cluster', i, '_GO.csv'))
  print(str(tmp))
  return(0)
})

pdf(paste0(pdir, '/hm_GO_term5.pdf'), width = 6.8, height = 3.5)
print(plotGOEnrich(goRes))
dev.off()

pdf(paste0(pdir, '/hm_GO_term10.pdf'), width = 6.8, height = 3.5)
print(plotGOEnrich(goRes, n = 10))
dev.off()

# --------------------------------------
# compare original and fitted expression
# --------------------------------------
# png(paste0(pdir, 'fitHm.png'),width = 4000,height = 2200,res = 300)
# plotFitHm(Res, type = 'variable')
# dev.off()
# 
# png(paste0(pdir, 'fitHm_rownames.png'),width = 12000,height = 10000,res = 300)
# print(plotFitHm(Res, type='variable', showRowName = T, cellWidthTotal = 1000, cellHeightTotal = length(Res$cluster) * 10))
# dev.off()


# png(paste0(pdir, 'DiffFitHm.png'),width = 4000,height = 2200,res = 250)
# plotDiffFitHm(Res, type = 'variable', cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE)
# dev.off()

# png(paste0(pdir, 'DiffFitHm3.png'),width = 4500,height = 2200,res = 250)
# plotDiffFitHm3(Res,  cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = FALSE, break.0 = FALSE)
# dev.off()

# png(paste0(pdir, 'DiffFitHm3_rownames.png'),width = 7500,height = 3500,res = 300)
# plotDiffFitHm3(Res, showRowName = T, cellWidthTotal = 300, cellHeightTotal = length(Res$cluster) * 10, sep = ':.*', subsampleCell = FALSE, break.0 = FALSE)
# dev.off()
# 
# png(paste0(pdir, 'DiffFitHm4_rownames.png'),width = 7500,height = 3500,res = 300)
# plotDiffFitHm4(Res, showRowName = T, cellWidthTotal = 300, cellHeightTotal = length(Res$cluster) * 10, sep = ':.*' , subsampleCell = FALSE, break.0 = FALSE)
# dev.off()

png(paste0(pdir, 'DiffFitHm5_rownames.png'),width = 12000,height = 3500,res = 300)
plotDiffFitHm5(Res, showRowName = T, cellWidthTotal = 300, cellHeightTotal = length(Res$cluster) * 10, sep = ':.*' , subsampleCell = FALSE, break.0 = FALSE)
dev.off()


groupdiff <- Res$covariateGroupDiff
mat1 <- groupdiff[names(clu)[clu == max(clu)-1], ]
mat2 <- groupdiff[names(clu)[clu == max(clu)], ]
dn1 <- dimnames(mat1)
dn2 <- dimnames(mat2)
v1 <- as.vector(mat1)
v2 <- as.vector(mat2)
mat1 <- matrix(v1/max(v1), nrow(mat1))
mat2 <- matrix(v2/max(v2), nrow(mat2))
dimnames(mat1) <- dn1
dimnames(mat2) <- dn2

d1 <- apply(mat1, 1, max) - apply(mat1, 1, min)
d2 <- apply(mat2, 1, max) - apply(mat2, 1, min)
mat1 <- mat1[order(d1), ]  
mat2 <- mat2[order(d2), ]  
Res$covariateGroupDiff <- rbind(groupdiff[names(clu)[clu %in% seq(1, max(clu)-2)], ], mat1, mat2)

png(paste0(pdir, 'DiffFitHm3_clu_type_changepoint_cor.png'),width = 5000,height = 3000,res = 100)
plotDiffFitHm3(Res, cellWidthTotal = 200, cellHeightTotal = 300, subsampleCell = F)
dev.off()

# ## ----------
# ## plot DDG
# ## ----------
# DDGType <- DDGType[diffgene]
# # id <- sort(sample(1:ncol(Res$populationFit[[1]]), ncol(Res$expr)))
# # Res$populationFit[[1]] <- Res$populationFit[[1]][, id]
# # Res$populationFit[[2]] <- Res$populationFit[[2]][, id]
# # colnames(Res$population[[1]]) <- colnames(Res$population[[2]]) <- id
# 
# for (i in unique(DDGType)){  ## debug -- ok!!
#   print(i)
#   gene <- names(DDGType)[DDGType == i]
#   png(paste0(pdir, 'diffgene_sampleFit_', i, '.png'), width = 4000, height = 2500, res = 200)
#   print(plotGene(Res, gene = gene[1:min(length(gene), 25)], plot.point = T, point.size = 0.1, variable = 'type'))
#   dev.off()
# }
# 
# for (i in unique(DDGType)){
#   print(i)
#   gene <- names(DDGType)[DDGType == i]
#   png(paste0(pdir, 'diffgene_groupFit_', i, '.png'), width = 2500, height = 2500, res = 200)
#   print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)], subSampleNumber=1000))
#   dev.off()
# }
# 
# for (i in 1:max(Res$cluster)){
#   print(i)
#   gene <- rownames(res)[res$cluster == i]
#   png(paste0(pdir, 'diffgene_groupFit_cluster', i, '.png'), width = 2500, height = 2500, res = 200)
#   print(plotGenePopulation(testobj = Res, type = 'variable', gene = gene[1:min(length(gene), 100)], subSampleNumber=1000))
#   dev.off()
# }

expr <- Res$expr[c(row.names(Res$expr)[1:100],"RPS4Y1:ENSG00000129824","CD3D:ENSG00000167286","CD3G:ENSG00000160654"),]
cellanno <- Res$cellanno  
pseudotime <- Res$pseudotime
design <- Res$design
permuiter=100;EMmaxiter=100;EMitercutoff=0.1;verbose=F;ncores=detectCores();test.type='Time';fit.resolution = 1000;return.all.data = TRUE;demean = FALSE;overall.only = F;test.method = 'permutation';ncores.fit = 1;fix.all.zero = TRUE;cutoff = 1e-5

rr <- testpt(Res$expr[c(row.names(Res$expr)[1:100],"RPS4Y1:ENSG00000129824","CD3D:ENSG00000167286","CD3G:ENSG00000160654"),], Res$cellanno, Res$pseudotime, Res$design,test.type='Variable',ncores=10)

png(paste0(pdir, 'RPS4Y1.png'),width = 1000,height = 1000,res = 100)
plotGene(Res, gene = "RPS4Y1:ENSG00000129824", plot.point = T)
dev.off()

png(paste0(pdir, 'RPS4Y1_new.png'),width = 1000,height = 1000,res = 100)
plotGene(rr, gene = "RPS4Y1:ENSG00000129824", plot.point = T)
dev.off()


