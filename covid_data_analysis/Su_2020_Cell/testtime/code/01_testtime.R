m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)

library(here)
setwd(here())
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')

meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_4levels.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]

source('function/01_function.R')
rdir <- paste0('covid/Su_2020_Cell/testtime/result/', m, '/')
system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design[,1,drop=F], test.type='Time', ncores = 16, demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), ncores.fit = 8)
})
saveRDS(res, paste0(rdir, 'testtime_res.rds'))


