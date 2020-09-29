## covid diffpt
library(parallel)
m <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/PRJCA002413_pbmc/data/proc/pt/expr.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/testpt.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/fitpt.R')
pt <- 1:ncol(m)
names(pt) <- colnames(m)
cellanno = data.frame(cell=colnames(m), sample = sub(':.*','',colnames(m)), stringsAsFactors = FALSE)
rownames(cellanno) <- cellanno[,1]
unis <- unique(cellanno[,2])
design <- cbind(1,as.numeric(grepl('Healthy',unis)))
rownames(design) <- unis
m <- m[rowMeans(m > 0.1) > 0.01,]
res <- testpt(expr=m, cellanno=cellanno, pseudotime=pt, design=design,type='Variable')

dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result', recursive = T)
saveRDS(res, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/testpt_res.rds')
