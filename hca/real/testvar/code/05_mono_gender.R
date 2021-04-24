m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)

library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/monocyte/'
rdir <- paste0('hca/real/testvar/result/', m, '/monocyte/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)
system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design,ncores=24, test.type = 'Variable', demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), cutoff = 1e-3)
})
saveRDS(res, paste0(rdir, 'gender_res.rds'))



