library(VGAM)
library(parallel)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/saverlog_pm.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/design.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/cellanno_pm.rds')

cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
x <- design[match(cellanno[,2],rownames(design)),2]

expressionFamily <- 'uninormal'
pval <- mclapply(1:nrow(d),function(i) {
  i <- d[i,]
  full_model_fit <- VGAM::vglm(i~x, epsilon=1e-1, family=expressionFamily)
  reduced_model_fit <- VGAM::vglm(i~1, epsilon=1e-1, family=expressionFamily) 
  
  lrt <- VGAM::lrtest(full_model_fit,reduced_model_fit) 
  pval=lrt@Body["Pr(>Chisq)"][2,]
},mc.cores=20)
pval <- unlist(pval)
names(pval) <- rownames(d)
fdr <- p.adjust(pval,method='fdr')
saveRDS(fdr,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/monocle2/fdr.rds')
