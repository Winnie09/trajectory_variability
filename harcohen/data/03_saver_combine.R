library(here)
setwd(here())
ddir <- 'harcohen/data/saverEachSample/'
af <- list.files(ddir)
f <- af[1]
data <- lapply(af, function(f){
  print(f)
  d <- readRDS(paste0(ddir, f))$estimate
  d <- log2(d + 1)
})
  
data <- do.call(cbind, data)
saveRDS(data, 'harcohen/data/saverEachSample/log2_saver.rds')



