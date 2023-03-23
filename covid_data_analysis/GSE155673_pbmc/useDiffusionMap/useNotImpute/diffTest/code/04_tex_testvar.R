library(here)
here()
cellanno <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','cellanno.rds'))
expr <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','log2norm.rds'))
pt <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','pseudotime.rds'))
ds <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','design.rds'))
design = data.frame(intercept = 1, type = ifelse(ds$type=='HD', 0, 1), stringsAsFactors = F)
rownames(design) <- rownames(ds)
pt <- pt[names(pt) %in% colnames(expr)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source(here('function','01_function.R'))
## trenddiff
rdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result')
dir.create(rdir, recursive = T)

system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 4, demean = FALSE)
})
saveRDS(res, paste0(rdir, '/tex_testvar.rds'))


