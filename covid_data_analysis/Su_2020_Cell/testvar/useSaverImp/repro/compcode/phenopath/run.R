library(phenopath)
seed <- as.numeric(commandArgs(trailingOnly = T)[1])
partid <- as.numeric(commandArgs(trailingOnly = T)[2])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

rdir <- paste0('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/phenopath10/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

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
expr <- expr[, cellanno[,1]]

set.seed(seed)
samp1 <- sample(rownames(design),nrow(design)/2)
samp2 <- setdiff(rownames(design),samp1)
cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
cellanno2 <- cellanno[cellanno[,2] %in% samp2,]
if (partid==1) {
  cellanno <- cellanno1
} else if (partid==2) {
  cellanno <- cellanno2
}

expr <- expr[,cellanno[,1]]

x <- cbind(design[match(cellanno[,2],rownames(design)),2],model.matrix(~cellanno[,2])[,-1])
print(dim(x))
fit <- phenopath(t(expr), x, elbo_tol = 1e-10, thin = 10,maxiter=10)
saveRDS(fit,file=paste0(rdir,seed,'_',partid,'_fit_res.rds'))
sig <- significant_interactions(fit)
saveRDS(sig,file=paste0(rdir,seed,'_',partid,'_sig_res.rds'))

