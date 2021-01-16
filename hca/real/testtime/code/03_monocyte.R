library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/monocyte/'
rdir <- 'hca/real/testtime/result/monocyte/'
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

m = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[,1,drop=FALSE]

system.time({
  res <- testpt(expr=m, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=1, permuiter=3, type = 'Time', demean = FALSE)
})
saveRDS(res, paste0(rdir, 'testtime_res.rds'))

