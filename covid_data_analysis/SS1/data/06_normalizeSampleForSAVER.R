library(here)
setwd(here())
rdir <- 'covid/SS1/data/eachSampleNorm/'
dir.create(rdir, showWarnings = F, recursive = T)
cnt <- readRDS('covid/SS1/data/count.rds')
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


