geneProp <- 0.2
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testvar/addMultiSignalUsingExpr/')
# ddir <- './testtime/data/data/'
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/null/'
rdir <-'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testvar/addMultiSignalUsingExpr/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')

### load saver, and count matrix
expr <- readRDS(paste0(ddir, 'hsc_mep_ery_saver.rds')) ## log2
cnt <- readRDS(paste0(ddir, 'hsc_mep_ery_count.rds'))

### prepare count, imputed, selected genes
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
rownames(expr) <- sapply(rownames(expr), function(i) sub('_','-',i))
expr.bak = expr
expr <- 2^expr - 1


design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

set.seed(12345)
selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))
othgene <- setdiff(rownames(expr), selgene)
saveRDS(selgene, paste0(rdir, 'selgene/selgene.rds'))

### order cells by pseudotime, and get the optimal clusters using saver imputed 80% genes
pt <- readRDS(paste0(ddir,'pseudotime.rds'))
expr <- expr[, pt[,1]] ##  NOT log2
tmp <- expr.bak[othgene, pt[,1]]  ##  log2
cnt <- cnt[, colnames(expr)]

### seperate by two group
sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
expr1 <- expr[, sample %in% rownames(design[design[,'group']==0, , drop = F])] ## NOT log2
expr2 <- expr[, sample %in% rownames(design[design[,'group']==1, , drop = F])]  ## NOT log2
tmp1 <- tmp[, colnames(expr1)] ## log2
tmp2 <- tmp[, colnames(expr2)] ## log2

### scale NOT log2 each sample 
library(matrixStats)
scalematrix <- function(data) {  ## standadize for rows
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  (data - cm) / csd
}

tmp1p <- sapply(colnames(tmp1), function(i) sub(':.*', '', i))
tmp1_list <- lapply(unique(tmp1p), function(i){
  scalematrix(tmp1[, tmp1p == i])
})
tmp1 <- do.call(cbind, tmp1_list)
tmp1 <- tmp1[, colnames(expr1)]
tmp1[is.na(tmp1)] <- 0

set.seed(12345)
clu <- kmeans(tmp1,10,iter.max = 1000)
clu <- clu$cluster
saveRDS(clu, paste0(rdir, 'clu/geneCluster.rds'))

# > table(clu)
# clu
#    1    2    3    4    5    6    7    8    9   10 
# 1726  437  890  580  652  861  600  189  823  498 

# set.seed(12345)
# clu2 <- mykmeans(t(tmp1))$cluster
# > table(clu2)
# clu2
#    1    2    3 
# 2991 1772 2493 

# -----
# plot
# -----
# # plot cluster mean of standadized expression
# clumean <- sapply(1:max(clu),function(i) colMeans(tmp1[clu==i,, drop=F]))
# clumean <- clumean[pt[,1][pt[,1] %in% colnames(tmp1)], ]
# rownames(clumean) <- 1:nrow(clumean)
# pd <- reshape2::melt(clumean)
# colnames(pd) <- c('pt','cluster','expr')
# 
# library(ggplot2)
# ggplot(pd,aes(x=pt,y=expr)) + 
#   geom_point(col='grey', size=0.01) + 
#   geom_smooth() + 
#   facet_wrap(~cluster,scales = 'free') + ylab('cluster mean of standadized expression') +
#   theme_classic()
# 
# # plot all genes expr vs. pt in each group with linew
# pd.clu = lapply(1:max(clu), function(i) reshape2::melt(tmp1[clu==i, ]))
# pd.clu = do.call(rbind, pd.clu)
# colnames(pd.clu) <- c('gene', 'cell', 'expr')
# pd.clu$clu <- clu[pd.clu$gene]
# pd.clu$pt <- pt[pd.clu$gene, 2]
# ggplot(data = pd.clu, aes(x = pt, y = expr, group = gene)) +
#   geom_line(color='grey', size = 0.01) +
#   theme_classic() +
#   facet_wrap(~clu, scales = 'free')
## plot
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


# -------------------------------------------------------
### add signal to permuted expression (both cnt and expr)
# -------------------------------------------------------
dir.create(paste0(rdir, 'count/'), showWarnings = F, recursive = T)
dir.create(paste0(rdir, 'saver/'), showWarnings = F, recursive = T)
for (j in seq(1,4)) { # signal from strongest -> weakest
  print(j)
  fromgene <- lapply(seq(1, max(clu)), function(i){
    print(i)
    clumat = expr1[names(clu[clu==i]), ]  ## NOT log2
    allp <- sapply(colnames(expr1), function(i) sub(':.*','', i))
    sds <- sapply(unique(allp), function(p){
      tmp <- clumat[, allp == p]
      sd <- apply(tmp, 1, sd) 
      rank(sd)
    })
    s <- names(sort(rowMeans(sds), decreasing = TRUE)) ## genes sd largest -> smallest
    
    set.seed(12345)
    fromgene <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], floor(length(selgene)/10), replace = T)
  })
  fromgene = unlist(fromgene)
  if (length(fromgene) < length(selgene)){
    fromgene = c(fromgene, fromgene[seq(1, length(selgene) - length(fromgene))])
  }
  saveRDS(fromgene, paste0(rdir, 'fromgene/', j, '.rds'))
  resexpr <-  expr1[selgene, ] + expr1[fromgene, ] # NOT log2
  mat <- rbind(resexpr[selgene, ], expr1[othgene,])
  mat <- cbind(mat, expr2[rownames(mat), ])
  saveRDS(mat, paste0(rdir, 'saver/', j, '.rds'))

  resexpr <-  cnt[selgene, colnames(expr1)] + cnt[fromgene, colnames(expr1)]
  mat <- rbind(resexpr[selgene, ], cnt[othgene, colnames(expr1)])
  mat <- cbind(mat, cnt[rownames(mat), colnames(expr2)])
  saveRDS(mat, paste0(rdir, 'count/', j, '.rds'))

}

  

    
    
