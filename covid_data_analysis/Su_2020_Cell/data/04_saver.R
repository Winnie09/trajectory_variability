library(SAVER)
library(Matrix)
library(parallel)
library(here)
packageVersion('SAVER')
setwd(here())
ddir <- 'covid/Su_2020_Cell/data/eachSampleNorm/'
rdir <- 'covid/Su_2020_Cell/data/eachSampleSaver/'
dir.create(rdir, recursive = TRUE)
f <- as.character(commandArgs(trailingOnly = T)[[1]])
print(f)
d <- readRDS(paste0(ddir, f))
d.saver <- saver(d, ncores = detectCores(), size.factor = 1)
saveRDS(d.saver, paste0(rdir, f))

