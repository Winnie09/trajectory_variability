library(here)
here()
setwd(here('covid','data','temra/'))
cellanno <- readRDS('cellanno.rds')
expr <- readRDS('log2norm.rds')
pt <- readRDS('pseudotime.rds')
ds <- readRDS('design.rds')
design = data.frame(intercept = 1, type = ifelse(ds$type=='HD', 0, 1), stringsAsFactors = F)
rownames(design) <- rownames(ds)
pt <- pt[names(pt) %in% colnames(expr)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source(here('function/01_function.R'))

rdir <- here('covid','result')
dir.create(rdir, recursive = T)
# testtime
# res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Time', ncores = 8, demean = TRUE)
# saveRDS(res, paste0(rdir, 'temra_testtime_centered.rds,'temra_testtime.rds'))

# testvar
## trenddiff
# res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 8, demean = TRUE)
# saveRDS(res, paste0(rdir, '/temra_trenddiff_centered,'/temra_trenddiff.rds'))

res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 8, demean = FALSE)
saveRDS(res, paste0(rdir,'/temra_NOT_centered', '/temra_NOT_centered.rds'))

# ## meandiff
# meanres <- meandiff(expr = expr, cellanno = cellanno, design = design, ncores =8)
# saveRDS(meanres, paste0(rdir, '/temra_meandiff','/temra_meandiff.rds'))
