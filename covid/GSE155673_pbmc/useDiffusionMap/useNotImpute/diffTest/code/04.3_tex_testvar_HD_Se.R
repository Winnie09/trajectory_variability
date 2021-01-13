library(here)
source(here('function','01_function.R'))
rdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', 'tex_testvar_HD_Se')
dir.create(rdir, recursive = T)

cellanno <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','cellanno.rds'))
expr <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','log2norm.rds'))
pt <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','pseudotime.rds'))
ds <- readRDS(here('covid','GSE155673_pbmc','data','useDiffusionMap','tex','design.rds'))
## select only HD and Se
design = data.frame(intercept = 1, type = ifelse(ds$type=='HD', 0, ifelse(ds$type == 'Se', 1, NA)), stringsAsFactors = F)
rownames(design) <- rownames(ds)
design <- design[complete.cases(design), ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
expr <- expr[, cellanno[,1]]
pt <- pt[names(pt) %in% colnames(expr)]

expr <- expr[rowMeans(expr>0.1)>0.01, ]
system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 4, demean = FALSE)
})
saveRDS(res, paste0(rdir, '/res.rds'))



