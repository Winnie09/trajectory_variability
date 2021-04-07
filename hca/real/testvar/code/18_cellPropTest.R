library(here)
setwd(here())
source('function/01_function.R')
for (path in c('lymph', 'erythroid', 'monocyte')){
  ddir <- paste0('hca/real/build_from_tree_variability/result/', path, '/')
  rdir <- paste0('hca/real/testvar/result/EM_pm/', path, '/')
  cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
  design = readRDS(paste0(ddir, 'input_design.rds'))
  pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))
  design = design[, 1:2]
  design[,2] <- ifelse(design[,2] == 'male', 0, 1)
  res <- cellPropTest(cellanno = cellanno, pseudotime = pseudotime, design = design, ncores = 4)
  saveRDS(res, paste0(rdir, 'cell_proportion_test_gender.rds'))
  
  design = readRDS(paste0(ddir, 'input_design.rds'))
  design = design[, c(1,3)]
  design[,2] <- as.numeric(design[,2])
  res <- cellPropTest(cellanno = cellanno, pseudotime = pseudotime, design = design, ncores = 4)
  saveRDS(res, paste0(rdir, 'cell_proportion_test_age.rds'))
}
