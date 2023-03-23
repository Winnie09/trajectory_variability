library(here)
setwd(here())
expr <- readRDS('covid/SS1/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/SS1/data/tActivate_pseudotime.rds')
cellanno <- readRDS('covid/SS1/data/cellanno.rds')
design <- readRDS('covid/SS1/data/design_numeric_Se_Mi.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]

source('function/01_function.R')
rdir <- paste0('covid/SS1/testvar/result/EM_pm/Se_Mi/')
dir.create(rdir, showWarnings = F, recursive = T)
system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, test.type='Variable', ncores = 25, demean = FALSE, test.method = 'permutation', ncores.fit = 8)
})
saveRDS(res, paste0(rdir, 'numeric_res.rds'))


