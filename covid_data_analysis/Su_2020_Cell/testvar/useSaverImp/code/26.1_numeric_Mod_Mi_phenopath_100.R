library(phenopath)
library(here)
setwd(here())
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/phenopath/Mod_Mi/')

expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
d <- expr[, cellanno[,1]]

####

d <- d[,names(pt)]
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
x <- cbind(design[match(cellanno[,2],rownames(design)),2],model.matrix(~cellanno[,2])[,-1])
print(dim(x))
fit <- phenopath(t(d), x, elbo_tol = 1e-10, thin = 100,maxiter=100)
saveRDS(fit,file=paste0(rdir,'fit_res.rds'))
sig <- significant_interactions(fit)
saveRDS(sig,file=paste0(rdir,'sig_res.rds'))

