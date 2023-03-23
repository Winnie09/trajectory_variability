m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)

library(here)
setwd(here())
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Se.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]

source('function/01_function.R')
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/', m, '/Mod_Se/')
dir.create(rdir, recursive = T, showWarnings = F)
system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, test.type='Variable', ncores = 12, demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'))
})
saveRDS(res, paste0(rdir, 'numeric_res.rds'))


