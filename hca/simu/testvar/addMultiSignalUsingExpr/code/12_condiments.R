library(tradeSeq)
f <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/count/')[as.numeric(commandArgs(trailingOnly = T))]
d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/count/',f))
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/null/pseudotime_pm_df.rds')
d <- d[,pt[,1]]
design = cbind(1,matrix(c(1,1,0,0,1,1,0,0), nrow=8))
rownames(design) = paste0('BM',seq(1,8))
cellanno = data.frame(cell=colnames(d), sample = sub(':.*','', colnames(d)), stringsAsFactors = FALSE)
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]

x <- design[match(cellanno[,2],rownames(design)),2]
x <- factor(x)
sce <- fitGAM(counts = d, pseudotime = matrix(1:ncol(d),ncol=1), cellWeights = matrix(1,nrow=ncol(d),ncol=1),condition=x, verbose = TRUE)
saveRDS(sce,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/condiments/sce_',f))
cond_genes <- conditionTest(sce)
saveRDS(cond_genes,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/condiments/cond_genes_',f))

