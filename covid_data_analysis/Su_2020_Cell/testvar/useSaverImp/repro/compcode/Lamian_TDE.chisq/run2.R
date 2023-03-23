seed <- as.numeric(commandArgs(trailingOnly = T))
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
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

source('function/01_function.R')
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/repro/compres/Lamian_TDE.chisq/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)
samp1 <- sample(rownames(design),nrow(design)/2)
samp2 <- setdiff(rownames(design),samp1)
cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
cellanno2 <- cellanno[cellanno[,2] %in% samp2,]

expr=expr[,cellanno2[,1]]
res2 <- testpt(expr=expr, cellanno=cellanno2, pseudotime=pt[cellanno2[,1]], design=design[samp2,1,drop=FALSE], test.type='Time', ncores = 15, demean = FALSE, test.method = 'chisq')
saveRDS(res2, paste0(rdir, seed,'_2.rds'))
