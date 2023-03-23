library(here)
setwd(here())
ddir <- 'hca/real/testvar/result/'
path = 'erythroid'
for (path in c('monocyte', 'lymph', 'erythroid')){
  print(path)
  res <- readRDS(paste0(ddir, path, '/age_fdr_res.rds'))
  source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
  rn <- sub(':.*', '', rownames(res))
  res <- res[!duplicated(rn), ]
  rownames(res) <- sub(':.*', '', rownames(res))
  print(dim(res[res[,1]<0.05, ]))
  print(rownames(res[res[,1]<0.05, ])[1:100])
  go <- myGO(rownames(res[res[,1]<0.05, ]), rownames(res))
  go <- go[order(go[,'FDR'], -abs(go[, 'FC'])), ]
  print(summary(go[, 'FDR']))
  print(head(go))
}

  
