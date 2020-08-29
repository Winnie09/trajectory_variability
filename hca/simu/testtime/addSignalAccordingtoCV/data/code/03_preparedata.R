geneProp <- 0.2

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testtime/data/data/'
### load saver, and count matrix
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
expr <- readRDS('./testtime/data/data/null/hsc_mep_ery_saver.rds')
cnt <- readRDS('./testtime/data/data/null/hsc_mep_ery_count.rds')

### prepare count, imputed, selected genes
expr <- 2^expr - 1
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
rownames(expr) <- sapply(rownames(expr), function(i) sub('_','-',i))
cnt <- cnt[, colnames(expr)]
allp = sub(':.*','', colnames(expr))
sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
set.seed(12345)
selgene <- sample(row.names(expr), round(geneProp * nrow(expr)))
othgene <- setdiff(rownames(expr), selgene)
saveRDS(selgene, paste0(rdir, 'selgene/selgene.rds'))

### permute the pseudotime for each sample
pmlist <- lapply(unique(allp), function(p){
  tmp <- expr[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmmat <- do.call(cbind, pmlist)

pmlist <- lapply(unique(allp), function(p){ ###### -->
  tmp <- cnt[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmcnt <- do.call(cbind, pmlist)
identical(colnames(pmmat), colnames(pmcnt)) ######### <------
identical(rownames(pmmat), rownames(pmcnt)) ######### <------

### order cells by pseudotime, and get the optimal clusters using saver imputed 80% genes
pt <- readRDS(paste0(rdir,'null/pseudotime.rds'))
expr <- expr[, pt[,1]]
tmp <- expr[othgene, pt[,1]]

library(matrixStats)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  (data - cm) / csd
}

set.seed(12345)
tmp <- scalematrix(tmp)
library(parallel)
# rss <- mclapply(1:20,function(clunum) {
#   tmp <- kmeans(tmp,clunum,iter.max = 1000)
#   tmp$betweenss/tmp$totss
# },mc.cores=20)
# rss <- unlist(rss)

clu <- kmeans(tmp,10,iter.max = 1000)
clu <- clu$cluster
clumean <- sapply(1:10,function(i) colMeans(tmp[clu==i,]))
clumean <- clumean[pt[,1],]
rownames(clumean) <- 1:nrow(clumean)
library(reshape2)
pd <- melt(clumean)
colnames(pd) <- c('pt','cluster','expr')
# ggplot(pd,aes(x=pt,y=expr)) + geom_point(col='grey', size=0.01) + geom_smooth() + facet_wrap(~cluster,scales = 'free') + ylab('scaled expression')
saveRDS(clu, paste0(rdir, 'null/geneCluster.rds'))


### add signal to permuted expression (both cnt and expr) 
dir.create(paste0(rdir, 'count/'), showWarnings = F, recursive = T)
dir.create(paste0(rdir, 'saver/'), showWarnings = F, recursive = T)

for( i in seq(1,10)){
  print(i)
  for (j in seq(1,4)) {
    # i = 9
    
    clumat = expr[names(clu[clu==i]), ]
    s <- names(sort(apply(clumat,1,sd)/rowMeans(clumat)))
    
    # j = 1:4
    set.seed(12345)
    addgene <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], length(selgene), replace = T)
    addexpr <- expr[addgene, pt[,1]]
    resexpr <-  pmmat[selgene, pt[,1]] + addexpr
    mat <- rbind(resexpr[selgene, pt[,1]], pmmat[othgene, pt[,1]])
    saveRDS(mat, paste0(rdir, 'saver/clusterType', i, '_', j,'.rds'))
    
    addexpr <- cnt[addgene, pt[,1]]
    resexpr <-  pmcnt[selgene, pt[,1]] + addexpr
    mat <- rbind(resexpr[selgene, pt[,1]], pmcnt[othgene, pt[,1]])
    saveRDS(mat, paste0(rdir, 'count/clusterType', i, '_', j,'.rds'))
  }
}
    

  
  ################################
  # # plot
  # par(mfrow=c(1,3))
  # plot(pt[,2],pmmat['HBB:ENSG00000244734', pt[,1]], cex=0.1) 
  # plot(pt[,2],expr['HBB:ENSG00000244734', pt[,1]], cex=0.1) 
  # plot(pt[,2],resexpr['HBB:ENSG00000244734', pt[,1]], col='red', cex=0.1)
  # 
  # par(mfrow=c(1,3))
  # plot(pt[,2],pmmat[selgene[2], pt[,1]], cex=0.1) 
  # plot(pt[,2],expr[addgene[2], pt[,1]], cex=0.1) 
  # plot(pt[,2],resexpr[selgene[2], pt[,1]], col='red', cex=0.1)
  # 
  # 
  # cv1 <- apply(resexpr[selgene, ], 1, sd)/rowMeans(resexpr[selgene, ])
  # cv2 <- apply(pmmat[othgene, ], 1, sd)/rowMeans(pmmat[othgene, ])
  # ggplot(data=rbind(data.frame(cv = cv1, g = 'selgene'), data.frame(cv = cv2, g = 'othgene')), aes(x=g, y = cv, col=g)) + geom_boxplot(outlier.colour = 'NA')
  # 
  # c1 <- sapply(selgene, function(g) cor(resexpr[g, ], seq(1, ncol(resexpr))))
  # c2 <- sapply(othgene, function(g) cor(pmmat[g, ], seq(1, ncol(pmmat))))
  # ggplot(data=rbind(data.frame(c = c1, g = 'selgene'), data.frame(c = c2, g = 'othgene')), aes(x=g, y = c, col=g)) + geom_boxplot(outlier.colour = 'NA')
  # 
  # 
  # mat <- rbind(resexpr[selgene, pt[,1]], pmmat[othgene, pt[,1]])
  # maxknot <- 5
  # knots = seq(1, ncol(mat),length.out=maxknot+2)[2:(maxknot+1)]
  # x <- cbind(1, bs(1:ncol(mat), knots = knots))
  # beta <- mat %*%  t(chol2inv(chol(crossprod(x))) %*% t(x))
  # fit <- t(x %*% t(beta))
  # ss1 <- sapply(rownames(mat), function(i){sum((fit[i, ] - mat[i, ])^2)})
  # ss2 <- sapply(rownames(mat), function(i){sum((mean(mat[i, ]) - mat[i, ])^2)})
  # 
  # f <- ((ss2-ss1)/ (9 - 1) )/(ss1/ (ncol(mat) - 9))
  # # full model rss/fulm dof
  # # dof = num.datapoint - num.prarmeters in full(9)
  # # null mdoel rss - full model ss / (fullmodel num.para - null model para)
  # ggplot() + geom_boxplot(data=rbind(data.frame(f = f[selgene], g = 'selgene'), data.frame(f = f[othgene], g = 'othgene')), aes(x=g, y = f, col=g), outlier.colour = 'NA') + ylim(c(-1,2))
  # 
  # 

