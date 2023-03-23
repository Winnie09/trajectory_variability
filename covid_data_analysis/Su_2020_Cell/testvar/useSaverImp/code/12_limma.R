library(limma)
library(here)
setwd(here())
source('function/01_function.R')
rdir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/limma/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

for (comparison in c('HD_Se', 'HD_Mi','Mod_Se','HD_Mod', 'Se_Mi', 'Recovered_Deceased')){
  print(comparison)
  expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
  pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
  meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
  cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')

  design <- readRDS(paste0('covid/Su_2020_Cell/data/design_numeric_', comparison, '.rds'))
  rownames(cellanno) <- cellanno[,1]
  cellanno <- cellanno[pt, ]
  cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
  pt <- pt[names(pt) %in% cellanno[,1]]
  expr <- expr[, names(pt)]
  expr <- expr[rowMeans(expr>0.1)>0.01, ]
  expr <- expr[, cellanno[,1]]
  agg <- vapply(unique(cellanno[,2]), function(i)
    rowMeans(expr[, cellanno[,2] == i, drop = FALSE], na.rm = TRUE), 
    numeric(nrow(expr)))
  agg <- agg[, rownames(design)]

  res <- topTable(eBayes(lmFit(agg, design)),n=nrow(agg),coef=2)
  res <- res[order(res[,5], -abs(res[,1])), ]
  saveRDS(res, paste0(rdir, comparison, '.rds'))
}

