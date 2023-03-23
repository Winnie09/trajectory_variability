m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)

library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/erythroid/'
rdir <- paste0('hca/real/testvar/result/', m, '/erythroid/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)
system.time({
  res <- testpt(expr=expr[expr=expr[sample(1:nrow(expr), 10), ], ], cellanno=cellanno, pseudotime=pseudotime, design=design,ncores=4, test.type = 'Variable', demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), EMitercutoff=0.1, permuiter = 4)
})
saveRDS(res, paste0(rdir, 'gender_res1e3.rds'))





