library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/lymph/'
rdir <- 'hca/real/testtime/result/lymph/'
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

m = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = matrix(1, nrow=length(unique(cellanno[,2])))
rownames(design) <- unique(cellanno[,2])
colnames(design) <- 'intercept'

system.time({
  res <- testpt(expr=m, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=4, type = 'Time', demean = FALSE)
})
saveRDS(res, paste0(rdir, 'testtime_res.rds'))

