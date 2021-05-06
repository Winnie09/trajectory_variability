library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/monocyte/'
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
rownames(cellanno) <- cellanno[,1]

design = readRDS(paste0(ddir, 'input_design.rds'))
design = design[,1,drop=FALSE]
design <- design[sort(rownames(design)), , drop=F]

for (id in 1:8){
  print(id)
  rdir <- paste0('hca/real/testtime_leave_one_out/data/data/leaveout_', id, '/')
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  selcell1 <- cellanno[cellanno[,2] %in% rownames(design)[id], 1]
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
  
  saveRDS(selcell1, paste0(rdir, 'selcell_1sample.rds'))
  saveRDS(design1, paste0(rdir, 'design_1sample.rds'))
  saveRDS(cellanno1, paste0(rdir, 'cellanno_1sample.rds'))
  saveRDS(pt1, paste0(rdir, 'pt_1sample.rds'))
  saveRDS(expr1, paste0(rdir, 'expr_1sample.rds'))
  
  saveRDS(selcell2, paste0(rdir, 'selcell_7sample.rds'))
  saveRDS(design2, paste0(rdir, 'design_7sample.rds'))
  saveRDS(cellanno2, paste0(rdir, 'cellanno_7sample.rds'))
  saveRDS(pt2, paste0(rdir, 'pt_7sample.rds'))
  saveRDS(expr2, paste0(rdir, 'expr_7sample.rds'))
}

