library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/erythroid/'
rdir <- 'hca/real/testvar/result/erythroid/'
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

m = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[, c(1,3)]
design[,2] <- as.numeric(design[,2])
system.time({
  res <- testpt(expr=m, cellanno=cellanno, pseudotime=pseudotime, design=design,ncores=1, permuiter=3, type = 'Variable', demean = FALSE)
})
saveRDS(res, paste0(rdir, 'age_res.rds'))

