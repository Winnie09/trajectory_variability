library(here)
setwd(here())
ddir <- 'covid/SS1/data/eachSampleSaver/'
rdir <- 'covid/SS1/data/'
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



