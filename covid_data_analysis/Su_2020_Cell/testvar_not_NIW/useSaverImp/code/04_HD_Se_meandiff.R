library(here)
setwd(here())
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
# design <- readRDS('covid/Su_2020_Cell/data/design_numeric_4levels.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_HD_Se.rds')
# design <- readRDS('covid/Su_2020_Cell/data/design_numeric_HD_others.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
str(expr)

source('function/01_function.R')
rdir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
expr <- expr[,cellanno[,1]]
agg <- vapply(unique(cellanno[,2]), function(i)
  rowMeans(expr[, cellanno[,2] == i, drop = FALSE], na.rm = TRUE), 
  numeric(nrow(expr)))
agg <- agg[, rownames(design)]
library(limma)
res <- topTable(eBayes(lmFit(agg, design)),n=nrow(agg),coef=2)
saveRDS(res, paste0(rdir, 'numeric_HD_Se_meandiff.rds'))

