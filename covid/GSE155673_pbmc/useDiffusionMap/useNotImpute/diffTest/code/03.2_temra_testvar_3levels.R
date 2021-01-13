library(here)
here()
cellanno <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','temra','cellanno.rds'))
expr <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','temra','log2norm.rds'))
pt <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','temra','pseudotime.rds'))
ds <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','temra','design.rds'))

design = data.frame(intercept = 1, type = ifelse(ds$type=='HD', 0, ifelse(ds$type == 'Mod', 1, 2)), stringsAsFactors = F)
rownames(design) <- rownames(ds)
pt <- pt[names(pt) %in% colnames(expr)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source(here('function','01_function.R'))
## trenddiff
rdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', 'temra_testvar_3levels')
dir.create(rdir, recursive = T)

system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 4, demean = FALSE)
})
saveRDS(res, paste0(rdir, '/res.rds'))


