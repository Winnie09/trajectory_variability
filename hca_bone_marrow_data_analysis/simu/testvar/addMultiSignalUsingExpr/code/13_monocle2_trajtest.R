library(VGAM)
library(parallel)
f <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/saver/')[as.numeric(commandArgs(trailingOnly = T))]
d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/saver/',f))
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/null/pseudotime_pm_df.rds')
d <- d[,pt[,1]]
design = cbind(1,matrix(c(1,1,0,0,1,1,0,0), nrow=8))
rownames(design) = paste0('BM',seq(1,8))
cellanno = data.frame(cell=colnames(d), sample = sub(':.*','', colnames(d)), stringsAsFactors = FALSE)
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
pt <- 1:ncol(d)
x <- design[match(cellanno[,2],rownames(design)),2]

pval <- mclapply(1:nrow(d),function(i) {
print(i)
  s <- d[i,]
  full_model_fit <- VGAM::vglm(s~x+sm.ns(pt, df=3), epsilon=1e-1, family='uninormal')
  reduced_model_fit <- VGAM::vglm(s~x, epsilon=1e-1, family='uninormal')
  lrt <- VGAM::lrtest(full_model_fit,reduced_model_fit) 
  lrt@Body["Pr(>Chisq)"][2,]
},mc.cores=20)
pval <- unlist(pval)
res <- data.frame(pval=pval,fdr=p.adjust(pval,method='fdr'))
rownames(res) <- rownames(d)
saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/monocle2_trajtest/',f))


