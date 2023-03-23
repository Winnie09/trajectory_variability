library(here)
setwd(here())
source('function/01_function.R')
for (comparison in c('Mod_Mi', 'HD_Mi', 'Mod_Se', 'HD_Se')){ ## 
  print(comparison)
  pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
  cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
  design <- readRDS(paste0('covid/Su_2020_Cell/data/design_numeric_', comparison, '.rds'))
  rownames(cellanno) <- cellanno[,1]
  cellanno <- cellanno[pt, ]
  cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
  pt <- pt[names(pt) %in% cellanno[,1]]
  
  rdir <- paste0('covid/Su_2020_Cell/testtime/result/EM_pm/cellPropTest/', comparison, '/')
  dir.create(rdir, recursive = T, showWarnings = F)
  system.time({
    res <- cellPropTest(cellanno=cellanno, pseudotime=pt, design=design, test.type = 'time')
  })
  saveRDS(res, paste0(rdir, 'cellPropTest_res.rds'))
}

