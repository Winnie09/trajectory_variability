library(here)
here()
cellanno <- readRDS(here('covid','data','tex','cellanno.rds'))
expr <- readRDS(here('covid','data','tex','log2norm.rds'))
pt <- readRDS(here('covid','data','tex','pseudotime.rds'))
ds <- readRDS(here('covid','data','tex','design.rds'))
design = data.frame(intercept = 1, type = ifelse(ds$type=='HD', 0, 1), stringsAsFactors = F)
rownames(design) <- rownames(ds)
pt <- pt[names(pt) %in% colnames(expr)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source(here('function/01_function.R'))

rdir <- here('covid','result')
dir.create(rdir, recursive = T)

# testtime
# res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Time', ncores = 8, demean = T)
# saveRDS(res, paste0(rdir, /tex_testtime_centered,'/tex_testtime.rds'))

# # trenddiff
# res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 8, demean = T)
# saveRDS(res, paste0(rdir, /tex_trenddiff_centered,'/tex_trenddiff.rds'))

res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 8, demean = FALSE)
saveRDS(res, paste0(rdir, '/tex_NOT_centered','/tex_NOT_centered.rds'))

# ## meandiff
# meanres <- meandiff(expr = expr, cellanno = cellanno, design = design, ncores =8)
# saveRDS(meanres, paste0(rdir, '/tex_meandiff,'/tex_meandiff.rds'))
