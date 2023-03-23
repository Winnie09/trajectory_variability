m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)

library(here)
setwd(here())
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/build_from_tree_variability/result/lymph/'
rdir <- paste0('hca/real/testtime/result/', m, '/lymph/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = matrix(1, nrow=length(unique(cellanno[,2])))
rownames(design) <- unique(cellanno[,2])
colnames(design) <- 'intercept'

system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=16, test.type = 'Time', demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), ncores.fit = 48)
})

saveRDS(res, paste0(rdir, 'testtime_res.rds'))

