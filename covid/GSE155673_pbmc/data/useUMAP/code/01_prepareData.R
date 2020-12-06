library(parallel)
library(Seurat)
pt = readRDS('/dcl02/hongkai/data/covid/temp/order.rds')
smeta <- read.table('/dcl02/hongkai/data/whou/covid/imp/data/meta/sampleMeta2.txt', stringsAsFactors = FALSE)
smeta[,1] <- sub('_', '', sub('_RNA', '', smeta[,1]))
rownames(smeta) <- smeta[,6]
colnames(smeta) <- c('samplename', 'type2','type3','age','gender','sample')
smeta$type <- sub('-.*', '', smeta[,2])
v1 <- sub(':.*', '', sub('GSE155673_', '', pt[[1]]))
v2 <- sub('.*:', '', sub('GSE155673_', '', pt[[1]]))
v1 <- smeta[match(v1,smeta[,1]),6]
pt[[1]] <- paste0(v1,':', v2)

expr <- sapply(unique(v1), function(i){
  print(i)
  tmp <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/imp/saver/res/sample/', i, '.rds'))$estimate
  tmp2 <- log2(tmp[,colnames(tmp) %in% pt[[1]]] + 1)
})
expr <- do.call(cbind, expr)  
meta <- data.frame(cell = colnames(expr), sample = sub(':.*', '', colnames(expr)), stringsAsFactors = FALSE)
saveRDS(meta, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/cellanno.rds')
saveRDS(expr, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/log2norm.rds')
pt1 = seq(1, length(pt[[1]]))
names(pt1) <- pt[[1]]
saveRDS(pt1, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/pseudotime.rds')
saveRDS(smeta, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/temra/design.rds')

v1 <- sub(':.*', '', sub('GSE155673_', '', pt[[2]]))
v2 <- sub('.*:', '', sub('GSE155673_', '', pt[[2]]))
v1 <- smeta[match(v1,smeta[,1]),6]
pt[[2]] <- paste0(v1,':', v2)
expr2 <- sapply(unique(v1), function(i){
  tmp <- readRDS(paste0('/dcl02/hongkai/data/whou/covid/imp/saver/res/sample/', i, '.rds'))$estimate
  tmp2 <- log2(tmp[,colnames(tmp) %in% pt[[2]]] + 1)
})
expr2 <- do.call(cbind, expr2)  
meta2 <- data.frame(cell = colnames(expr2), sample = sub(':.*', '', colnames(expr2)), stringsAsFactors = FALSE)
saveRDS(meta2, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/tex/cellanno.rds')
saveRDS(expr2, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/tex/log2norm.rds')
pt2 <- seq(1, length(pt[[2]]))
names(pt2) <- pt[[2]]
saveRDS(pt2, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/tex/pseudotime.rds')
saveRDS(smeta, '/dcl02/hongkai/data/whou/trajectory_variability/covid/data/tex/design.rds')



