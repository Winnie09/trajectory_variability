library(phenopath)
f <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/saver/')[as.numeric(commandArgs(trailingOnly = T))]
d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/saver/',f))
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/null/pseudotime_pm_df.rds')
d <- d[,pt[,1]]
design = cbind(1,matrix(c(1,1,0,0,1,1,0,0), nrow=8))
rownames(design) = paste0('BM',seq(1,8))
cellanno = data.frame(cell=colnames(d), sample = sub(':.*','', colnames(d)), stringsAsFactors = FALSE)

cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
x <- cbind(design[match(cellanno[,2],rownames(design)),2],model.matrix(~cellanno[,2])[,-1])
print(dim(x))
fit <- phenopath(t(d), x, elbo_tol = 1e-10, thin = 500,maxiter=500)
saveRDS(fit,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/phenopath500/fit_',f))
sig <- significant_interactions(fit)
saveRDS(sig,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/phenopath500/sig_',f))

