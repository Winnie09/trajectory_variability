rm(list=ls())
geneProp <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
addSignalType <- as.character(commandArgs(trailingOnly = TRUE)[[2]])
addSignalPara <-  as.numeric(commandArgs(trailingOnly = TRUE)[[3]])
## be reminded that the pca is done on cv > 0.5 (for all samples) genes. You should redo it !!!
# geneProp <- 0.1
# addSignalType <- 'linear'
# addSignalPara <-  5

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
library(parallel)
library(splines)
library(limma)

## get expr
expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
allp = sub('_.*','', colnames(expr))
# cv <- apply(expr,1,sd)/rowMeans(expr)
# allcv <- sapply(unique(allp), function(i){
#   tmp <- expr[, grepl(i, colnames(expr))]
#   cv <- apply(tmp,1,sd)/rowMeans(tmp)
#   tmp <- tmp[cv>0.5,]
#   rownames(tmp)
# })
# # allcv <- sapply(unique(allp), function(i){
# #   tmp <- expr[, grepl(i, colnames(expr))]
# #   findVariableGene(tmp)
# # })
# allcv <- unique(unlist(allcv))
# allcv <- allcv[!is.na(allcv)]
# expr <- expr[allcv, ]
# saveRDS(expr,'./addsignal/result/expr.rds')
# expr = readRDS('./addsignal/result/expr.rds')


## get pca
##dr <- prcomp(t(expr),scale=T)$x
##dr <- dr[,1:2]
# dr <- prcomp(t(expr[cv>0.5,]),scale=T)$x
# dr <- dr[,1:2]
# saveRDS(dr, './addsignal/result/dr_trendVar.rds') 
dr <- readRDS('./addsignal/result/dr.rds')

# library(ggplot2)
# ggplot(data.frame(x=dr[,1],y=dr[,2],ct=ct[id]),aes(x=x,y=y,col=ct)) + geom_point() + theme_classic()

## change cells'  sample label
# v1 = sapply(colnames(expr), function(i){
#   sub(':.*','',i)
# })
# v2 = sapply(colnames(expr), function(i){
#   paste0(strsplit(i,':')[[1]][2:3], collapse=':')
# })
# set.seed(12345)
# v = paste0(sample(v1), ':', v2)
# colnames(expr) <- rownames(dr) <- v

library(TSCAN)
pseudotime <- TSCANorder(exprmclust(t(dr),reduce = F,clusternum=4,clustermethod='kmeans'),orderonly = T)
# saveRDS(pseudotime,'./addsignal/result/pseudotime.rds')
psn <- 1:length(pseudotime)
names(psn) <- pseudotime
order = data.frame(Pseudotime = seq(1, length(pseudotime)), Cell=pseudotime)


sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

## one group along pseudotime
rtgene <- sapply(c(1,0), function(i){
  print(i)
  selectsample <- rownames(design[design[,1]==i,,drop=F])
  print(str(selectsample))
  selectcell <- names(sample[sample %in% selectsample])
  res <- OneGroupDiffALongTime(expr = expr[, selectcell], sample=sample[selectcell], dr = dr[selectcell,], pseudotime=pseudotime[pseudotime %in% selectcell], num.base = 3) 
  rownames(res[res[,3]<0.05,])
})

expr <- expr[intersect(rtgene[[1]], rtgene[[2]]), ] ##  [1:2795, 1:13269]  ##  [1:6234, 1:13269]

set.seed(12345)
selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))
# r2 <- exprdiff(expr=expr,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE,fstatonly=FALSE,permutime=1000)
# library(pROC)
# true = sapply(rownames(expr), function(i) ifelse(i %in% selgene, 1, 0))
# auc(roc(true, r[,'P.Value']))

# set.seed(12345)
# expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=0)
# r <- exprdiff(expr=expr,design=design,sample=sample, dr=dr,pseudotime=pseudotime,permutation=FALSE, num.base=2)
# r <- r[order(r[,3]),]
# dim(r[r[,3]<0.05, ])
# siggene <- rownames(r[r[,3]<0.05, ])

# dir.create(paste0('./addsignal/result/limma/',addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
# df <- t(sapply(0:20, function(i){
#   print(i)
#   expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=i)
#   r <- exprdiff(expr=expr2,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=FALSE)
#   saveRDS(r, paste0('./addsignal/result/limma/',addSignalType,'/',geneProp,'_',i,'.rds'))
#   r <- r[order(r[,3]),]
#   sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
#   res <- list()
#   res[['res']] <- r
#   res[['sensfdr']] <- c(i,AreaUnderSensFdr(sensfdr))
#   print(c(i,AreaUnderSensFdr(sensfdr)))
#   return(c(i,AreaUnderSensFdr(sensfdr)))
# }))
# 
# colnames(df)[1] <- c('Parameter')
# saveRDS(df, paste0('./addsignal/result/limma/',addSignalType,'/',geneProp,'_Para_FdrDiff_Area_Group1256.rds'))


dir.create(paste0('./addsignal/result/permu/',addSignalType,'/'), showWarnings = FALSE, recursive = TRUE)
expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=addSignalPara)
res <- exprdiff(expr=expr2,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE, permutime=10000)
r <- res[['res']]
r <- r[order(r[,3]),]
sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
res[['sensfdr']] <- c(i,AreaUnderSensFdr(sensfdr))
saveRDS(res, paste0('./addsignal/result/permu/', addSignalType,'/', geneProp,'_',addSignalPara,'.rds'))



# ## plot top differential genes 
# library(gridExtra)
# library(ggplot2)
# expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = addSignalType, parameter=3)
# r <- exprdiff(expr=expr2,design=design,sample=sample, dr=dr,pseudotime=pseudotime,permutation=FALSE)
# r <- r[order(r[,3]),]
# plist <- plotGeneGenderDiff(GeneSet = rownames(r)[1:16], mat = exprs, order = order, sample=sample, design=design, statistics = r)
# grid.arrange(grobs=plist,nrow=4)

# plot one gene 
# g <- 'ARHGAP23:ENSG00000275832'
# plist <- plotGeneGenderDiff(GeneSet = g, mat = expr2, order = order, sample=sample, design=design, statistics = r)
# grid.arrange(grobs=plist,nrow=1)
# 
#     linedlist <- lapply(unique(pd$Patient)[1:7], function(p){
#       print(p)
#       tmat = mat[g,grepl(p,colnames(mat)),drop=F]
#       trainX = order[match(colnames(tmat), order$Cell),'Pseudotime']      ### use time 
#       pred <- get_spline_fit(tmat, trainX=trainX, fit.min=min(order$Pseudotime), fit.max=max(order$Pseudotime))
#       tmpdf <- data.frame(Expr=pred[1,], Pseudotime=trainX, Patient=p, Group = ifelse(design[names(unlist(sapply(rownames(design), function(i) grep(i, p)))),1]==1,GroupNames[1],GroupNames[2]))
      
# df <- t(sapply(1:10, function(i){
#   print(i)
#   expr2 <- AddSignal(expr = expr, sample = sample, SelectGene = selgene, SelectSample = paste0('BM',c(1,2,5,6)), pseudotime = pseudotime, method = 'linear', parameter=i)
#   r <- exprdiff(expr=expr2,design=design,sample=sample,dr=dr,pseudotime=pseudotime,permutation=TRUE, permutime=10000)
#   saveRDS(r, paste0('./addsignal/result/permu/linear/',i,'.rds'))
#   r <- r[['res']]
#   r <- r[order(r[,3]),]
#   sensfdr <- SensFdr(Order = rownames(r), TruePositive = selgene, statistics=r)
#   print(c(i,AreaUnderSensFdr(sensfdr)))
#   c(i,AreaUnderSensFdr(sensfdr))
# }))
# colnames(df)[1] <- c('parameter')
# saveRDS(df, './addsignal/result/Group1256_Linear_Para_FdrDiff_Area_permutation.rds')
 


rm(list=ls())


