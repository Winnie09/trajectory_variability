m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)

library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/monocyte/'
rdir <- paste0('hca/real/testtime/result/', m, '/monocyte/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[,1,drop=FALSE]

system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=24, test.type = 'Time', demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), ncores.fit = 48)
})
saveRDS(res, paste0(rdir, 'testtime_res.rds'))


