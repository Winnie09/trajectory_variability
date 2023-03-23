# --------------
# global setting
# --------------
method <- as.character(commandArgs(trailingOnly = T)[[1]])
print(method)

suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
library(here)
setwd(here())
ddir <- 'hca/real/build3traj/manual/result/lymph/'
rdir <- paste0('hca/real/testtime/result/', method, '/lymph/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')
expr = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pseudotime = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = matrix(1, nrow=length(unique(cellanno[,2])))
rownames(design) <- unique(cellanno[,2])
colnames(design) <- 'intercept'

# ----------------------------------
# test  [one group along pseudotime]
# ----------------------------------

if (method == 'tscan'){
  expr <- expr[, names(pseudotime)]
  res <- TSCAN_time(expr=expr,pseudotime=pseudotime)
}

if (method == 'monocle2'){
  expr <- expr[, names(pseudotime)]
  res <- monocle2_time(expr=expr,pseudotime=pseudotime)
}

saveRDS(res, paste0(rdir, 'testtime_res.rds'))

