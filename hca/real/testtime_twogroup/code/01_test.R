library(here)
setwd(here())
for (path in c('monocyte', 'lymph', 'erythroid')){
  ddir <- paste0('hca/real/build_from_tree_variability/result/', path, '/')
  rdir <- paste0('hca/real/testtime_twogroup/data/', path, '/')
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
  design1 <- design[c(1,2,5,6),,drop=F]
  design2 <- design[setdiff(rownames(design), rownames(design1)),,drop=F]
  cellanno1 <- cellanno[selcell1, ]
  cellanno2 <- cellanno[selcell2, ]
  pseudotime <- pseudotime[names(sort(pseudotime))]
  pt1 <- seq(1, length(selcell1))
  names(pt1) <- selcell1
  pt2 <- seq(1, length(selcell2))
  names(pt2) <- selcell2
  
  saveRDS(selcell1, paste0(rdir, 'selcell1.rds'))
  saveRDS(design1, paste0(rdir, 'design1.rds'))
  saveRDS(cellanno1, paste0(rdir, 'cellanno1.rds'))
  saveRDS(pt1, paste0(rdir, 'pseudotime1.rds'))
  saveRDS(expr[, selcell1], paste0(rdir, 'expr1.rds'))
  
  saveRDS(selcell2, paste0(rdir, 'selcell2.rds'))
  saveRDS(design2, paste0(rdir, 'design2.rds'))
  saveRDS(cellanno2, paste0(rdir, 'cellanno2.rds'))
  saveRDS(pt2, paste0(rdir, 'pseudotime2.rds'))
  saveRDS(expr[,selcell2], paste0(rdir, 'expr2.rds'))
}

