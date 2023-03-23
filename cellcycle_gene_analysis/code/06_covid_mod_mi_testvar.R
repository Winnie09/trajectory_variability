library(here)
setwd(here())
d  = readRDS('cellcycle/res/covid_Mod_Mi_cellcycle_score.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds') ## Mod = 0, Mi = 1

## subset mod and mi
cellanno = cellanno[cellanno[,2] %in% rownames(design), ]
selcell = cellanno[, 1]
# mat = mat[, selcell]
pt = pt[selcell]

source('function/01_function.R')
res = testpt(expr = t(d[, 1:3]), design = design, cellanno = cellanno, pseudotime = pt, test.type = 'Variable', ncores = 1)

saveRDS(res, 'cellcycle/res/covid_Mod_Mi_testres.rds')


