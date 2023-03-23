library(here)
setwd(here())
expr <- readRDS('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/data/matrix/saver.rds') ## [1:22401, 1:32819] 
for (path in c('lymph', 'monocyte', 'erythroid')){
  rdir <- paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/real/traj/result/', path, '/')
  print(rdir)
  pt <- readRDS(paste0(rdir, 'pseudotime.rds'))
  m <- expr[, pt]
  m <- m[rowMeans(m > 0.1) > 0.01, ]
  cellanno = data.frame(cell = colnames(m), sample = sub('_.*','', colnames(m)), stringsAsFactors = FALSE)
  design = data.frame(intercept = 1,          
                      gender = sub('.*:','', unique(cellanno$sample)),
                      age = sapply(unique(cellanno$sample), function(i) strsplit(i, ':')[[1]][2]),
                      stringsAsFactors = FALSE)
  rownames(design) <- unique(cellanno$sample)
  print(design)
  pseudotime = seq(1, length(pt))
  names(pseudotime) <- pt
  saveRDS(m, paste0(rdir,'input_expr.rds'))
  saveRDS(cellanno, paste0(rdir,'input_cellanno.rds'))
  saveRDS(design, paste0(rdir,'input_design.rds'))
  saveRDS(pseudotime, paste0(rdir,'input_pseudotime.rds'))
}

