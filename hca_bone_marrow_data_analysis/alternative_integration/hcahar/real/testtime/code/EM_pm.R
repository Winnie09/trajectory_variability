library(here)
setwd(here())
for (ct in list.files('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/real/traj/result/')) {
  ddir <- paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/real/traj/result/',ct)
  rdir <- paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/real/testtime/result/',ct)
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  source('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/function/01_function.R')
  
  m = readRDS(paste0(ddir, '/input_expr.rds'))
  cellanno = readRDS(paste0(ddir, '/input_cellanno.rds'))
  design = readRDS(paste0(ddir, '/input_design.rds'))
  pseudotime = readRDS(paste0(ddir, '/input_pseudotime.rds'))
  
  design = matrix(1, nrow=length(unique(cellanno[,2])))
  rownames(design) <- unique(cellanno[,2])
  colnames(design) <- 'intercept'
  
  system.time({
    res <- testpt(expr=m, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=20, test.type = 'Time', demean = FALSE)
  })
  saveRDS(res, paste0(rdir, 'EM_pm.rds'))
}

