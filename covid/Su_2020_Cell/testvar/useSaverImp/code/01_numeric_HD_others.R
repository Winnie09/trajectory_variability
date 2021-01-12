library(here)
setwd(here())
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_HD_others.rds') ##
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]

source('function/01_function.R')
rdir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
system.time({
  res <- testpt.seed(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 8, demean = FALSE)
})
saveRDS(res, paste0(rdir, 'numeric_HD_others_res.rds'))

