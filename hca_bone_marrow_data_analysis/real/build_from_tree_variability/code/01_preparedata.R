library(here)
setwd(here())
res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tree_variability/auto_pc_auto_nclu_module/result/infer_tree_structure_res.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/build_from_tree_variability/result/'  
expr <- readRDS('./hca/data/proc/matrix/saver.rds') ## [1:22401, 1:32819] 

for (path in c('lymph', 'monocyte', 'erythroid', 'clp')){
  print(path)
  rdir2 <- paste0(rdir, path, '/')
  dir.create(rdir2, showWarnings = F, recursive = T)
  if (path == 'lymph') {
    ord <- res$order[[1]]
  } else if (path == 'monocyte'){
    ord <- res$order[[2]]
  } else if (path == 'lymph'){
    ord <- res$order[4]
  } else if (path == 'clp'){
    ord <- res$order[[3]]
  }
  pt <- seq(1, length(ord))
  names(pt) <- ord
  m <- expr[, ord]
  m <- m[rowMeans(m > 0.1) > 0.01, ]
  cellanno = data.frame(cell = colnames(m), sample = sub('_.*','', colnames(m)), stringsAsFactors = FALSE)
  design = data.frame(intercept = 1,          
                      gender = sub('.*:','', unique(cellanno$sample)),
                      age = sapply(unique(cellanno$sample), function(i) strsplit(i, ':')[[1]][2]),
                      stringsAsFactors = FALSE)
  rownames(design) <- unique(cellanno$sample)
  print(design)
  
  saveRDS(m, paste0(rdir2,'input_expr.rds'))
  saveRDS(cellanno, paste0(rdir2,'input_cellanno.rds'))
  saveRDS(design, paste0(rdir2,'input_design.rds'))
  saveRDS(pt, paste0(rdir2,'input_pseudotime.rds'))
}


