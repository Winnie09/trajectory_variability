library(Matrix)
library(parallel)
library(here)
setwd(here())
rdir <- 'covid/Su_2020_Cell/data/eachSampleNorm/'
dir.create(rdir, recursive = T)
cnt <- readRDS(paste0(rdir, 'count.rds'))
ap <- sub(':.*', '', colnames(cnt))
libsize <- colSums(cnt)
libsize <- libsize/median(libsize) ####

nn <- sapply(unique(ap), function(p){
  print(p)
  tmp <- cnt[, ap == p, drop = FALSE]
  tmp <- sweep(tmp, 2, libsize[ap == p], '/')
  saveRDS(tmp, paste0(rdir, p, '.rds'))
  return(0)
})

