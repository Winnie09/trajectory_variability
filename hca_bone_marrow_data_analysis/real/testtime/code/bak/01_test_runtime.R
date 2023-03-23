m <- as.character(commandArgs(trailingOnly = T)[[1]])
print(m)

library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/lymph/'
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

m = 'EM_pm'
system.time({
  res <- testpt(expr=expr[1:3, ], cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=16, test.type = 'Time', demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), ncores.fit = 24)
})
#    user  system elapsed 
# 130.555  19.917  19.070 

m = 'EM_chisq'
system.time({
  res <- testpt(expr=expr[1:3, ], cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=16, test.type = 'Time', demean = FALSE, test.method = ifelse(m == 'EM_pm', 'permutation', 'chisq'), ncores.fit = 24)
})
  #  user  system elapsed 
  # 2.222   0.939   2.541 


