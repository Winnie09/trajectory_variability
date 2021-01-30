library(here)
setwd(here('hca','data','simu','testtime','addMultiSignalUsingExpr'))
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
dir.create('fromgene', recursive = TRUE, showWarnings = FALSE)

### load saver, count matrix, and pseudotime
geneProp <- 0.2
saverlog <- readRDS(here('hca','data','simu','testtime','poolSampleSignal','null','hsc_mep_ery_saver.rds'))
cnt <- readRDS(here('hca','data','simu','testtime','poolSampleSignal','null','hsc_mep_ery_count.rds'))
pt <- readRDS(here('hca','data','simu','testtime','poolSampleSignal','null','pseudotime.rds'))
saverlog <- saverlog[, pt[,1]]
cnt <- cnt[, pt[,1]]

### prepare count, imputed, selected genes
savercnt <- 2^saverlog - 1
rownames(savercnt) <- sapply(rownames(savercnt), function(i) sub('_','-',i))
rownames(saverlog) <- rownames(savercnt)
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
allp = sub(':.*','', colnames(savercnt))
names(allp) <- colnames(savercnt)
set.seed(12345)
selgene <- sample(row.names(savercnt), round(geneProp * nrow(savercnt)))
othgene <- setdiff(rownames(savercnt), selgene)
dir.create('selgene/', showWarnings = F, recursive = T)
saveRDS(selgene, './selgene/selgene.rds')

### permute the pseudotime for each sample
pmlist <- lapply(unique(allp), function(p){
  tmp <- savercnt[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmsavercnt <- do.call(cbind, pmlist)
pmsavercnt <- pmsavercnt[, pt[,1]]

pmlist <- lapply(unique(allp), function(p){ 
  tmp <- cnt[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmcnt <- do.call(cbind, pmlist)
pmcnt <- pmcnt[, pt[,1]]
identical(colnames(pmsavercnt), colnames(pmcnt)) 
identical(rownames(pmsavercnt), rownames(pmcnt)) 

## select highly pseudotime-variable genes 
library(splines)
tmp <- saverlog[othgene, ]
xm <- bs(1:ncol(tmp))
fstat <- t(sapply(rownames(tmp), function(i) { ## larger f, stronger signal
  summary(lm(tmp[i,pt[,1]]~xm))$fstatistic
}))

tmp.fit <- get_spline_fit(trainData = tmp, trainX = seq(1, ncol(tmp)), fit.min = 1, fit.max = ncol(tmp), fit.num.points = 1000, num.base=10, remove.correlated=TRUE)

cor(apply(tmp.fit, 1, sd),rowMeans(tmp.fit))
tmp.fit.sd <- apply(tmp.fit, 1, sd)
tmp.fit.mean <- rowMeans(tmp.fit)

xm <- bs(tmp.fit.mean)
resid <- resid(lm(tmp.fit.sd ~ xm))
pdf('./null/bs_fitted_sd_versus_mean_marked_by_residule.pdf', width = 3.8, height = 3.6)
plot(apply(tmp.fit, 1, sd)~rowMeans(tmp.fit), pch = 20, col = ifelse(resid>0, 'red', 'black'), xlab='mean', ylab='sd')
dev.off()
pdf('./null/bs_fitted_sd_versus_fstat_marked_by_residule.pdf', width = 3.8, height = 3.6)
plot(apply(tmp.fit, 1, sd)~fstat[,1], pch = 20, col = ifelse(resid>0, 'red', 'black'), main = paste0('PCC=',round(cor(apply(tmp.fit, 1, sd),fstat[,1]),2)), xlab='f statistics', ylab='sd')
dev.off()

othgene <- names(resid[resid >0])
saveRDS(othgene, './selgene/othgene.rds')

### get the optimal clusters using log2-saver imputed 80% genes
library(matrixStats)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  (data - cm) / csd
}

tmp <- saverlog[othgene, ]
tmp <- scalematrix(tmp)
library(parallel)
rss <- mclapply(1:10,function(clunum) {
  tmp <- kmeans(tmp,clunum,iter.max = 1000)
  tmp$betweenss/tmp$totss
},mc.cores=10)
rss <- unlist(rss)
names(rss) <- paste0(seq(1,10), 'cluster')
saveRDS(rss, './null/number_clusters_rss.rds')

pdf('./null/number_of_clusters_rss.pdf', width = 3.8, height = 3.8)
plot(rss,pch=20,xlab='number of clusters')
dev.off()

set.seed(12345)
clu <- kmeans(tmp,5,iter.max = 1000) ### 
clu <- clu$cluster
dir.create('null/', showWarnings = F, recursive = T)
saveRDS(clu, './null/geneCluster.rds')

clumean <- sapply(1:max(clu),function(i) colMeans(tmp[clu==i,]))
clumean <- clumean[pt[,1],]
rownames(clumean) <- 1:nrow(clumean)
library(ggplot2)
library(reshape2)
pd <- melt(clumean)
colnames(pd) <- c('pt','cluster','saverlog')
pdf('./null/gene_cluster.pdf', width = 6, height = 5)
ggplot(pd,aes(x=pt,y=saverlog)) + geom_point(col='grey', size=0.01) + geom_smooth() + facet_wrap(~cluster,scales = 'free') + ylab('Scaled gene expression') + xlab('Pseudotime') + theme_classic()
dev.off()

### add signal to permuted expression (both cnt and savercnt) 
dir.create('count/', showWarnings = F, recursive = T)
dir.create('saver/', showWarnings = F, recursive = T)
for (j in seq(1,4)) { # signal from 1 weakest to 4 strongest
  print(j)
  fromgene <- lapply(seq(1, max(clu)), function(i){
    print(i)
    clumat = saverlog[names(clu[clu==i]), ]  # log2 expression
    fstat.clu <- fstat[rownames(clumat), 1]
    s <- names(sort(fstat.clu))
    set.seed(12345)
    fromgenetmp <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], floor(length(selgene)/10), replace = T)
  }) ## j = 1, smallest fstat, weakest signal
  fromgene = unlist(fromgene)
  if (length(fromgene) < length(selgene)){
    set.seed(12345)
    fromgene = c(fromgene, sample(fromgene, length(selgene) - length(fromgene), replace = TRUE))
  }
  saveRDS(fromgene, paste0('fromgene/', j, '.rds'))
  
  resexpr <-  pmsavercnt[selgene, pt[,1]] + savercnt[fromgene, pt[,1]] # NOT log2
  mat <- rbind(resexpr[selgene, ], pmsavercnt[setdiff(rownames(savercnt),selgene),pt[,1]])
  saveRDS(mat, paste0('saver/', j,'.rds'))
  
  resexpr <-  pmcnt[selgene, pt[,1]] + cnt[fromgene, pt[,1]]
  mat <- rbind(resexpr[selgene, ], pmcnt[setdiff(rownames(cnt),selgene), pt[,1]])
  saveRDS(mat, paste0('count/', j,'.rds'))
  
  # ## double check
  # logmat <- log2(mat+1)
  # summary(sapply(sample(selgene,10),function(i) {  ## should be larger
  #   summary(lm(logmat[i,pt[,1]]~I(1:ncol(logmat))))$fstatistic[1]
  # }))
  # summary(sapply(sample(othgene,10),function(i) { ## should be much smaller
  #   summary(lm(logmat[i,pt[,1]]~I(1:ncol(logmat))))$fstatistic[1]
  # }))
}





