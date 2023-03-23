library(phenopath)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/case/phenopath10/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/case_design.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/expr.rds')

x <- cbind(design[match(cellanno[,2],rownames(design)),2],model.matrix(~cellanno[,2])[,-1])
print(dim(x))
fit <- phenopath(t(expr), x, elbo_tol = 1e-10, thin = 3,maxiter=3)
saveRDS(fit,file=paste0(rdir,'fit_res.rds'))
sig <- significant_interactions(fit)
saveRDS(sig,file=paste0(rdir,'sig_res.rds'))


