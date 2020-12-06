setwd('/dcl02/hongkai/data/whou/trajectory_variability/covid/data/tex/')
cellanno <- readRDS('cellanno.rds')
expr <- readRDS('log2norm.rds')
pt <- readRDS('pseudotime.rds')
ds <- readRDS('design.rds')
design = data.frame(intercept = 1, type = ifelse(ds$type=='HD', 0, 1), stringsAsFactors = F)
rownames(design) <- rownames(ds)
pt <- pt[names(pt) %in% colnames(expr)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source('/dcl02/hongkai/data/whou/trajectory_variability/function/01_function.R')

rdir <- '/dcl02/hongkai/data/whou/trajectory_variability/covid/result/'
dir.create(rdir, recursive = T)
# testtime
res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Time', ncores = 8)
saveRDS(res, paste0(rdir, 'tex_testtime.rds'))

