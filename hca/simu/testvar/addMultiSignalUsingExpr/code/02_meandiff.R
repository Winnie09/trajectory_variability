setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca')
rdir <- './simu/testvar/addMultiSignalUsingExpr/result/'
ddir <- './simu/testvar/addMultiSignalUsingExpr/data/'
method = 'meandiff'

for (signal in 1:4){
  print(signal)
  fn <- paste0(rdir, method,'/', signal, '.rds')
  print(fn)
  if (file.exists(fn)) break 
  suppressMessages(library(limma))
  source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
  
  ### prepare data
  expr <- readRDS(paste0(ddir, 'saver/', signal, '.rds'))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr>0)>0.01, ]
  cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
  design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
  
  agg <- vapply(unique(cellanno[,2]), function(i)
    rowMeans(expr[, cellanno[,2] == i, drop = FALSE], na.rm = TRUE), 
    numeric(nrow(expr)))
  agg <- agg[, paste0('BM', 1:8)]
  res <- topTable(eBayes(lmFit(agg, design)),n=nrow(agg),coef=2)
  saveRDS(res, fn)
}
  

