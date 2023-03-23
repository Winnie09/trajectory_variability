library(here)
packageVersion('SAVER')
setwd(here())
ddir <- 'harcohen/data/'
rdir <- 'harcohen/data/eachSample/'
dir.create(rdir, recursive = TRUE)
mat <- readRDS(paste0(ddir, 'norm.rds'))
expr <- 2^mat - 1
ap <- sub(':.*', '', colnames(mat))
nn <- sapply(unique(ap), function(p){
  print(p)
  tmp <- expr[, ap == p, drop = FALSE]
  saveRDS(tmp, paste0(rdir, p, '.rds'))
  return(0)
})
