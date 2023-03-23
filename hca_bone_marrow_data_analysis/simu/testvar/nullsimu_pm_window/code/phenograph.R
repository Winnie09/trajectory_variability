library(phenopath)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/saverlog_pm.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/design.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/cellanno_pm.rds')

cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
x <- cbind(design[match(cellanno[,2],rownames(design)),2],model.matrix(~cellanno[,2])[,-1])
print(dim(x))
fit <- phenopath(t(d), x, elbo_tol = 3e-6, thin = 40)
saveRDS(fit,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/phenopath/fit.rds')
sig <- significant_interactions(fit)
saveRDS(sig,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/phenopath/sig.rds')


