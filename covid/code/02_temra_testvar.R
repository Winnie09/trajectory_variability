cellanno <- readRDS('/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/cellanno.rds')
expr <- readRDS('/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/log2norm.rds')
pt <- readRDS('/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/pseudotime.rds')
ds <- readRDS('/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/design.rds')
design = data.frame(intercept = 1, type = ds$type, stringsAsFactors = F)
rownames(design) <- rownames(ds)
pt <- pt[names(pt) %in% colnames(expr)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source('/dcl02/hongkai/data/whou/trajectory_variability/function/01_function.R')
## trenddiff
rdir <- '/dcl02/hongkai/data/whou/trajectory_variability/covid/result/'
dir.create(rdir, recursive = T)


# res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Time', ncores = 8)
# saveRDS(res, paste0(rdir, 'temra_testtime.rds'))

res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 8)
saveRDS(res, paste0(rdir, 'temra_trenddiff.rds'))

# ------------------------
## meandiff
meanres <- meandiff(expr = expr, cellanno = cellanno, design = design, ncores =8)
saveRDS(meanres, paste0(rdir, 'temra_meandiff.rds'))

