library(here)
setwd(here())
source('./function/01_function.R')
for (path in c('lymph/', 'erythroid/', 'monocyte/')){
  print(path)
  ddir <- paste0('hca/real/build3traj/manual/result/', path)
  rdir <- paste0('hca/real/testvar/result/limma/', path)
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  print(rdir)
  
  m = readRDS(paste0(ddir, 'input_expr.rds'))
  cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
  design.bak = readRDS(paste0(ddir, 'input_design.rds'))
  pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
  
  design = design.bak[, 1:2]
  design[,2] <- ifelse(design[,2] == 'male', 0, 1)
  print(design)  
  res <- meandiff(expr = m, cellanno = cellanno, design = design)
  res <- res[order(res[,'adj.P.Val'], -abs(res[, 'logFC'])), ]
  print(sum(res[,'adj.P.Val']<0,05))
  saveRDS(res, paste0(rdir, 'gender_res.rds'))
  
  design = design.bak[, c(1,3)]
  design[,2] <- as.numeric(design[,2])
  print(design)
  res <- meandiff(expr = m, cellanno = cellanno, design = design)
  res <- res[order(res[,'adj.P.Val'], -abs(res[, 'logFC'])), ]
  print(sum(res[,'adj.P.Val']<0,05))
  saveRDS(res, paste0(rdir, 'age_res.rds'))
}
  



