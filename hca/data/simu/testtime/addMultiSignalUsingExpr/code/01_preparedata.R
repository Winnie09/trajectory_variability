## 03_preparedata.R
library(here)
setwd(here('hca','data','simu','testtime','addMultiSignalUsingExpr'))
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
dir.create('fromgene', recursive = TRUE, showWarnings = FALSE)

### load saver, and count matrix
geneProp <- 0.2
saverlog <- readRDS(here('hca','data','simu','testtime','poolSampleSignal','null','hsc_mep_ery_saver.rds'))
cnt <- readRDS(here('hca','data','simu','testtime','poolSampleSignal','null','hsc_mep_ery_count.rds'))

### prepare count, imputed, selected genes
savercnt <- 2^saverlog - 1
rownames(savercnt) <- sapply(rownames(savercnt), function(i) sub('_','-',i))
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
cnt <- cnt[, colnames(savercnt)]
allp = sub(':.*','', colnames(savercnt))
sample <- sub(':.*','',colnames(savercnt))
names(sample) <- colnames(savercnt)
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

pmlist <- lapply(unique(allp), function(p){ ###### -->
  tmp <- cnt[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmcnt <- do.call(cbind, pmlist)
identical(colnames(pmsavercnt), colnames(pmcnt)) ######### <------
identical(rownames(pmsavercnt), rownames(pmcnt)) ######### <------

### order cells by pseudotime, and get the optimal clusters using saver imputed 80% genes
pt <- readRDS(here('hca','data','simu','testtime','poolSampleSignal','null','pseudotime.rds'))
savercnt <- savercnt[, pt[,1]]
cnt <- cnt[, colnames(savercnt)]
allp <- sapply(colnames(savercnt), function(i) sub(':.*','', i))
tmp <- savercnt[othgene, pt[,1]]

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

# clumean <- sapply(1:10,function(i) colMeans(tmp[clu==i,]))
# clumean <- clumean[pt[,1],]
# rownames(clumean) <- 1:nrow(clumean)
# library(reshape2)
# pd <- melt(clumean)
# colnames(pd) <- c('pt','cluster','savercnt')
# ggplot(pd,aes(x=pt,y=savercnt)) + geom_point(col='grey', size=0.01) + geom_smooth() + facet_wrap(~cluster,scales = 'free') + ylab('scaled saver count expression')


### add signal to permuted expression (both cnt and savercnt) 
dir.create('count/', showWarnings = F, recursive = T)
dir.create('saver/', showWarnings = F, recursive = T)
for (j in seq(1,4)) { # signal from strongest -> weakest
  print(j)
  fromgene <- lapply(seq(1, max(clu)), function(i){
    print(i)
    clumat = savercnt[names(clu[clu==i]), ]  ## NOT log2

    sds <- sapply(unique(allp), function(p){
      tmp <- clumat[, allp == p]
      sd <- apply(tmp, 1, sd) 
      rank(sd)
    })
    s <- names(sort(rowMeans(sds), decreasing = TRUE)) ## genes sd largest -> smallest
    
    set.seed(12345)
    fromgenetmp <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], floor(length(selgene)/10), replace = T)
  }) ## j = 1, sd largest genes
  fromgene = unlist(fromgene)
  if (length(fromgene) < length(selgene)){
    fromgene = c(fromgene, fromgene[seq(1, length(selgene) - length(fromgene))])
  }
  saveRDS(fromgene, paste0('fromgene/', j, '.rds'))
  
  resexpr <-  pmsavercnt[selgene, ] + savercnt[fromgene, ] # NOT log2
  mat <- rbind(resexpr[selgene, ], savercnt[othgene,])
  saveRDS(mat, paste0('saver/', j,'.rds'))

  resexpr <-  pmcnt[selgene, colnames(savercnt)] + cnt[fromgene, colnames(savercnt)]
  mat <- rbind(resexpr[selgene, ], cnt[othgene, ])
  saveRDS(mat, paste0('count/', j,'.rds'))
}


