geneProp <- 0.2
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
ddir <- './testtime/data/data/'
rdir = './testvar/data/data/'
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')

### load saver, and count matrix
expr <- readRDS('./testtime/data/data/null/hsc_mep_ery_saver.rds') ## log2
cnt <- readRDS('./testtime/data/data/null/hsc_mep_ery_count.rds')

### prepare count, imputed, selected genes
expr <- 2^expr - 1
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
rownames(expr) <- sapply(rownames(expr), function(i) sub('_','-',i))
expr.bak = expr
cnt <- cnt[, colnames(expr)]
# allp = sub(':.*','', colnames(expr))
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

set.seed(12345)
selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))
othgene <- setdiff(rownames(expr), selgene)
saveRDS(selgene, paste0(rdir, 'selgene/selgene.rds'))

### order cells by pseudotime, and get the optimal clusters using saver imputed 80% genes
pt <- readRDS(paste0(ddir,'null/pseudotime.rds'))
expr <- expr[, pt[,1]] ## NOT log2
tmp <- expr.bak[othgene, pt[,1]]  ## log2


### seperate by two group
sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
expr1 <- expr[, sample %in% rownames(design[design[,'group']==0, , drop = F])]
expr2 <- expr[, sample %in% rownames(design[design[,'group']==1, , drop = F])]
tmp1 <- tmp[, colnames(expr1)]
tmp2 <- tmp[, colnames(expr2)]

### 
library(matrixStats)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  (data - cm) / csd
}


# tmp1 <- scalematrix(tmp1)
tmp1p <- sapply(colnames(tmp1), function(i) sub(':.*', '', i))
tmp1_list <- lapply(unique(tmp1p), function(i){
  scalematrix(tmp1[, tmp1p == i])
})
tmp1 <- do.call(cbind, tmp1_list)
tmp1 <- tmp1[, colnames(expr1)]
tmp1[is.na(tmp1)] <- 0
# library(parallel)
set.seed(12345)
clu <- kmeans(tmp1,10,iter.max = 1000)
clu <- clu$cluster

clumean <- sapply(1:10,function(i) colMeans(tmp1[clu==i,, drop=F]))
clumean <- clumean[pt[,1][pt[,1] %in% colnames(tmp1)], ]
rownames(clumean) <- 1:nrow(clumean)
library(reshape2)
pd <- melt(clumean)
colnames(pd) <- c('pt','cluster','expr')
#ggplot(pd,aes(x=pt,y=expr)) + geom_point(col='grey', size=0.01) + geom_smooth() + facet_wrap(~cluster,scales = 'free') + ylab('scaled expression')
# clumean <- sapply(1:10,function(i) colMeans(expr1[names(clu[clu==i]),, drop=F]))
# clumean <- clumean[pt[,1][pt[,1] %in% colnames(tmp1)], ]
# # rownames(clumean) <- 1:nrow(clumean)
# rownames(clumean) <- colnames(expr1) ##
# pd <- melt(clumean)
# v = sapply(as.character(pd[,1]), function(i) sub(':.*', '', i))
# rownames(clumean) <- 1:nrow(clumean)
# pd = melt(clumean)
# pd = cbind(pd, sample = v)
# colnames(pd) <- c('pt','cluster','expr', 'sample')
# ggplot(pd,aes(x=pt,y=expr,col=sample)) + geom_point(size=0.01) + geom_smooth() + facet_wrap(~cluster,scales = 'free') + ylab('scaled expression')
saveRDS(clu, paste0(rdir, 'null/geneCluster.rds'))

### add signal to permuted expression (both cnt and expr) 
dir.create(paste0(rdir, 'count/'), showWarnings = F, recursive = T)
dir.create(paste0(rdir, 'saver/'), showWarnings = F, recursive = T)
for( i in seq(1,10)){
  print(i)
  for (j in seq(1,4)) {
    clumat = expr1[names(clu[clu==i]), ]
    # clup <- sapply(colnames(clumat), function(i) sub(':.*', '', i))
    # g = rownames(clumat)[2]
    # ggplot(data = data.frame(x = seq(1, ncol(clumat)), y = clumat[g, ], clup = clup), aes(x = x, y = y, color=clup)) + geom_point( size = 0.1) 
    # 
    allp <- sapply(colnames(expr1), function(i) sub(':.*','', i))
    # s <- names(sort(apply(clumat,1,sd)/rowMeans(clumat))) ####
    sds <- sapply(unique(allp), function(p){
      tmp <- clumat[, allp == p]
      sd <- apply(tmp, 1, sd) 
      rank(sd)
    })
    s <- names(sort(rowMeans(sds), decreasing = TRUE))
    
    # j = 1:4
    set.seed(12345)
    addgene <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], length(selgene), replace = T)
    addexpr <- expr1[addgene, ]
    resexpr <-  expr1[selgene, ] + addexpr
    mat <- rbind(resexpr[selgene, ], expr1[othgene,])
    mat <- cbind(mat, expr2[rownames(mat), ])
    saveRDS(mat, paste0(rdir, 'saver/clusterType', i, '_', j,'.rds'))
    saveRDS(list(selgene = rownames(resexpr), addgene = rownames(addexpr)), paste0(rdir, 'saver/clusterType', i, '_', j,'_selgene_addgene.rds'))
    
    addexpr <- cnt[addgene, colnames(expr1)]
    resexpr <-  cnt[selgene, colnames(expr1)] + addexpr
    mat <- rbind(resexpr[selgene, ], cnt[othgene, colnames(expr1)])
    mat <- cbind(mat, cnt[rownames(mat), colnames(expr2)])
    saveRDS(mat, paste0(rdir, 'count/clusterType', i, '_', j,'.rds'))
  }
}


