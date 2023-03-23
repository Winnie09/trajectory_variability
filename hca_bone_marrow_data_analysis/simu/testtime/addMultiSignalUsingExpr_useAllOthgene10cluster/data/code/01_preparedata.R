## 03_preparedata.R
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
# sample <- sub(':.*','',colnames(savercnt))
# names(sample) <- colnames(savercnt)

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

pmlist <- lapply(unique(allp), function(p){ ###### -->
  tmp <- cnt[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmcnt <- do.call(cbind, pmlist)
pmcnt <- pmcnt[, pt[,1]]
identical(colnames(pmsavercnt), colnames(pmcnt)) ######### <------
identical(rownames(pmsavercnt), rownames(pmcnt)) ######### <------

### get the optimal clusters using log2-saver imputed 80% genes
tmp <- saverlog[othgene, ]
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
dir.create('null/', showWarnings = F, recursive = T)
saveRDS(clu, './null/geneCluster.rds')

clumean <- sapply(1:10,function(i) colMeans(tmp[clu==i,]))
clumean <- clumean[pt[,1],]
rownames(clumean) <- 1:nrow(clumean)
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
    fstat <- sapply(rownames(clumat), function(i) { ## larger f, stronger signal
      summary(lm(clumat[i,pt[,1]]~I(1:ncol(clumat))))$fstatistic[1]
  })
    names(fstat) <- rownames(clumat)
    s <- names(sort(fstat))
    set.seed(12345)
    fromgenetmp <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], floor(length(selgene)/10), replace = T)
  }) ## j = 1, smallest fstat, weakest signal
  fromgene = unlist(fromgene)
  if (length(fromgene) < length(selgene)){
    fromgene = c(fromgene, fromgene[seq(1, length(selgene) - length(fromgene))])
  }
  saveRDS(fromgene, paste0('fromgene/', j, '.rds'))
  
  resexpr <-  pmsavercnt[selgene, pt[,1]] + savercnt[fromgene, pt[,1]] # NOT log2
  mat <- rbind(resexpr[selgene, ], pmsavercnt[othgene,pt[,1]])
  saveRDS(mat, paste0('saver/', j,'.rds'))

  resexpr <-  pmcnt[selgene, pt[,1]] + cnt[fromgene, pt[,1]]
  mat <- rbind(resexpr[selgene, ], pmcnt[othgene, pt[,1]])
  saveRDS(mat, paste0('count/', j,'.rds'))
  
  # ## double check
  # logmat <- log2(mat+1)
  # summary(sapply(sample(selgene,10),function(i) {  ## should be larger
  #   summary(lm(logmat[i,pt[,1]]~I(1:ncol(logmat))))$fstatistic[1]
  # }))
  # summary(sapply(sample(othgene,10),function(i) { ## should be much smaller
  #   summary(lm(logmat[i,pt[,1]]~I(1:ncol(logmat))))$fstatistic[1]
  # }))
  # ##
}

