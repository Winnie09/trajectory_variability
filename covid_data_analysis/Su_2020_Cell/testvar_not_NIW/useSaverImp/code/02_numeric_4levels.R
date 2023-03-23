library(here)
setwd(here())
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_4levels.rds') ##
# design <- readRDS('covid/Su_2020_Cell/data/design_numeric_HD_Se.rds')
# design <- readRDS('covid/Su_2020_Cell/data/design_numeric_HD_others.rds')
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]
source('function/01_function.R')
rdir <- 'covid/Su_2020_Cell/testvar/useSaverImp/result/'
i = commandArgs(trailingOnly = T)[[1]][1]
if (i == '1'){
  system.time({
    res <- testpt(expr=expr[1:4000, ], cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 20, demean = FALSE)
  })
  saveRDS(res, paste0(rdir, 'numeric_4levels_res_1.rds'))
} else {
  system.time({
    res <- testpt(expr=expr[4001:nrow(expr), ], cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 20, demean = FALSE)
  })
  saveRDS(res, paste0(rdir, 'numeric_4levels_res_2.rds'))
}


  


