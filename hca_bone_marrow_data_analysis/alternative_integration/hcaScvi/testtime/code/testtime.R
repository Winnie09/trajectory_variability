library(here)
setwd(here())
oid <- as.numeric(commandArgs(trailingOnly = T))
expr <- readRDS('./hca/data/proc/matrix/saver.rds') ## [1:22401, 1:32819] 
ord <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/data/order.rds')
on <- names(ord)[oid]
ord <- ord[[oid]]
expr <- expr[,ord]
cellanno = data.frame(cell = colnames(expr), sample = sub('_.*','', colnames(expr)), stringsAsFactors = FALSE)
pseudotime = seq(1, length(ord))
names(pseudotime) <- ord

expr <- expr[rowMeans(expr > 0.1) > 0.1,]
source('./function/01_function.R')

design = matrix(1, nrow=length(unique(cellanno[,2])))
rownames(design) <- unique(cellanno[,2])
colnames(design) <- 'intercept'

system.time({
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, ncores=10, test.type = 'Time', demean = FALSE)
})
saveRDS(res, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/testtime/res/',on,'.rds'))



