library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/monocyte/'
rdir <- 'hca/real/testvar/result/monocyte/'
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

m = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)
system.time({
  res <- testpt(expr=m, cellanno=cellanno, pseudotime=pseudotime, design=design,ncores=12, type = 'Variable', demean = FALSE)
})
saveRDS(res, paste0(rdir, 'gender_res.rds'))

