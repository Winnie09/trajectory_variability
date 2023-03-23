library(phenopath)
library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/lymph/'
rdir <- paste0('hca/real/testvar/result/phenopath100/lymph/gender/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

d = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pt = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)
####

d <- d[,names(pt)]
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
x <- cbind(design[match(cellanno[,2],rownames(design)),2],model.matrix(~cellanno[,2])[,-1])
print(dim(x))
fit <- phenopath(t(d), x, elbo_tol = 3e-6, thin = 100,maxiter=100)
saveRDS(fit,file=paste0(rdir,'fit_res.rds'))
sig <- significant_interactions(fit)
saveRDS(sig,file=paste0(rdir,'sig_res.rds'))


