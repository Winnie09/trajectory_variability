library(here)
setwd(here())
ddir <- 'covid/Su_2020_Cell/data/eachSampleSaver/'
rdir <- 'covid/Su_2020_Cell/data/'
af <- list.files(ddir)
f <- af[1]
res <- lapply(af, function(f){
  print(paste0(which(af == f), '/', length(af)))
  m <- readRDS(paste0(ddir, f))
  m <- m$estimate
  m <- log2(m + 1)
})
  
mat <- do.call(cbind, res)
saveRDS(mat, paste0(rdir, 'saver_log2norm.rds'))
