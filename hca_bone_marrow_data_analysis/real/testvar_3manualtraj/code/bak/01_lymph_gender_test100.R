m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)
i = as.numeric(commandArgs(trailingOnly = T)[[2]])
print(paste0('The 100 truck: ', i))

library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/lymph/'
rdir <- paste0('hca/real/testvar/result/', m, '/lymph/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)
selectgene <- rownames(expr)[seq( (i-1)*100+1, min(i*100, nrow(expr) ))]
print(selectgene)

system.time({
  res <- testpt(expr=expr[selectgene, ], cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=12, test.type = 'Variable', test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), demean = FALSE)
})
saveRDS(res, paste0(rdir, 'gender_res', i, '.rds'))



