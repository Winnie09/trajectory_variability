library(here)
setwd(here())
expr <- readRDS('covid/Su_2020_Cell/data/log2norm.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_4levels.rds')
# design <- readRDS('covid/Su_2020_Cell/data/design_numeric_HD_Se.rds')
# design <- readRDS('covid/Su_2020_Cell/data/design_numeric_HD_others.rds')
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source('function/01_function.R')
rdir <- 'covid/Su_2020_Cell/testvar/result/'
system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 4, demean = FALSE)
})
saveRDS(res, paste0(rdir, 'numeric_4levels_res.rds'))
