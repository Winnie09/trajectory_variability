library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/monocyte/'
rdir <- 'hca/real/testtime_twogroup/data/data/'
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
rownames(cellanno) <- cellanno[,1]
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
design = design[,1,drop=FALSE]
design <- design[sort(rownames(design)), , drop=F]

selcell1 <- cellanno[cellanno[,2] %in% rownames(design)[c(1,2,5,6)], 1]
selcell2 <- setdiff(cellanno[,1], selcell1)
design1 <- design[c(1,2,5,6),, drop = F]
design2 <- design[setdiff(rownames(design), rownames(design1)),, drop = F]
cellanno1 <- cellanno[selcell1, ]
cellanno2 <- cellanno[selcell2, ]
pseudotime <- pseudotime[names(sort(pseudotime))]
pt1 <- pseudotime[selcell1]
pt2 <- pseudotime[selcell2]
expr1 <- expr[, selcell1]
expr2 <- expr[, selcell2]


saveRDS(selcell1, paste0(rdir, 'selcell1.rds'))
saveRDS(design1, paste0(rdir, 'design1.rds'))
saveRDS(cellanno1, paste0(rdir, 'cellanno1.rds'))
saveRDS(pt1, paste0(rdir, 'pt1.rds'))
saveRDS(expr1, paste0(rdir, 'expr1.rds'))

saveRDS(selcell2, paste0(rdir, 'selcell2.rds'))
saveRDS(design2, paste0(rdir, 'design2.rds'))
saveRDS(cellanno2, paste0(rdir, 'cellanno2.rds'))
saveRDS(pt2, paste0(rdir, 'pt2.rds'))
saveRDS(expr2, paste0(rdir, 'expr2.rds'))


