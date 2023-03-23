library(here)
setwd(here('hca','simu','testvar','addMultiSignalUsingExpr','data'))
source(here('function','01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
dir.create('fromgene', recursive = TRUE, showWarnings = FALSE)

### load saver, count matrix, and pseudotime
geneProp <- 0.2
saverlog <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','hsc_mep_ery_saver.rds'))
cnt <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','hsc_mep_ery_count.rds'))
pt <- readRDS(here('hca','simu','testtime','poolSampleSignal','data','null','pseudotime.rds'))
saverlog <- saverlog[, pt[,1]]
cnt <- cnt[, pt[,1]]

### permute the sample-cell relationsihp 
identical(colnames(saverlog), colnames(cnt))
sn <- sub(':.*', '', colnames(cnt))
cn <- sapply(1:ncol(cnt), function(i) sub(paste0(sn[i], ':'), '', colnames(cnt)[i]))
# plot(saverlog[grep('XIST', rownames(cnt)),], col = as.factor(sn %in% c(paste0('BM', c(1,2,5,6)))))
# plot(saverlog[grep('XIST', rownames(cnt)),], col = as.factor(sn))
set.seed(1234567)
pt[,1] <- colnames(saverlog) <- colnames(cnt) <- paste0(sample(sn),':', cn)
pt.v <- pt[,2]
names(pt.v) <- pt[,1]
saveRDS(pt, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/null/pseudotime_pm_df.rds')
saveRDS(pt.v, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/null/pseudotime_pm.rds')

### prepare count, imputed, selected genes
savercnt <- 2^saverlog - 1
rownames(savercnt) <- sapply(rownames(savercnt), function(i) sub('_','-',i))
rownames(saverlog) <- rownames(savercnt)
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
allp = sub(':.*','', colnames(savercnt))
names(allp) <- colnames(savercnt)
# BM1:52:female   BM2:50:male   BM3:39:male   BM4:29:male   BM5:29:male 
#          1453           513           794          2123          4085 
# BM6:26:female BM7:36:female BM8:32:female 
#          1633          1090          1578 
## do not change the order of group1 and group2 codes
savercnt2 <- savercnt[ , allp %in% paste0('BM', c(3,4,7,8))] ## group2, do not add signals
saverlog2 <- saverlog[ , allp %in% paste0('BM', c(3,4,7,8))] ## group2, do not add signals
cnt2 <- cnt[ , allp %in% paste0('BM', c(3,4,7,8))] ## group2, do not add signals
savercnt <- savercnt[ , allp %in% paste0('BM', c(1,2,5,6))] ## group1, add signals
saverlog <- saverlog[ , allp %in% paste0('BM', c(1,2,5,6))] ## group1, add signals
cnt <- cnt[ , allp %in% paste0('BM', c(1,2,5,6))] ## group1, add signals
allp = sub(':.*','', colnames(savercnt))
names(allp) <- colnames(savercnt)
pt1 <- pt[pt[,1] %in% colnames(cnt), ] ## pseudotime for only group1

set.seed(1234567)
selgene <- sample(row.names(savercnt), round(geneProp * nrow(savercnt)))
othgene <- setdiff(rownames(savercnt), selgene)
dir.create('selgene/', showWarnings = F, recursive = T)
saveRDS(selgene, './selgene/selgene.rds')

selgene1 <- selgene[1:round(length(selgene)/3)]
selgene2 <- selgene[(length(selgene1)+1) : round(2*length(selgene)/3+1)]
selgene3 <- selgene[(length(selgene1)+length(selgene2) + 1) : length(selgene)]
saveRDS(selgene1, './selgene/selgene1.rds')
saveRDS(selgene2, './selgene/selgene2.rds')
saveRDS(selgene3, './selgene/selgene3.rds')

## select highly pseudotime-variable genes 
library(splines)
tmp <- saverlog[othgene, ]
xm <- bs(1:ncol(tmp))
fstat <- t(sapply(rownames(tmp), function(i) { ## larger f, stronger signal
  summary(lm(tmp[i,pt1[,1]]~xm))$fstatistic
  # summary(lm(tmp[i,pt1[,1]]~xm[pt1[,1] %in% colnames(tmp),]))$fstatistic
  # summary(lm(tmp[i,pt1[,1]]~pt[which(pt1[,1] %in% pt[,1]),2]))$fstatistic
}))

tmp.fit <- get_spline_fit(trainData = tmp, trainX = seq(1, ncol(tmp)), fit.min = 1, fit.max = ncol(tmp), fit.num.points = 1000, num.base=10, remove.correlated=TRUE)

cor(apply(tmp.fit, 1, sd),rowMeans(tmp.fit))  ##  0.5540323
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
plot(rss,pch=20,xlab='number of clusters', ylab = 'betweenss/totss')
dev.off()

set.seed(1234567)
clu <- kmeans(tmp,5,iter.max = 1000) ### 
clu <- clu$cluster
dir.create('null/', showWarnings = F, recursive = T)
saveRDS(clu, './null/geneCluster.rds')

clumean <- sapply(1:max(clu),function(i) colMeans(tmp[clu==i,]))
clumean <- clumean[pt1[,1],]
rownames(clumean) <- 1:nrow(clumean)
library(ggplot2)
library(reshape2)
pd <- melt(clumean)
colnames(pd) <- c('pt','cluster','saverlog')
pdf('./null/gene_cluster.pdf', width = 6, height = 5)
ggplot(pd,aes(x=pt,y=saverlog)) + geom_point(col='grey', size=0.01) + geom_smooth() + facet_wrap(~cluster,scales = 'free') + ylab('Scaled gene expression') + xlab('Pseudotime') + theme_classic()
dev.off()

### permuted a group1 matrix for adding mean
allp <-  sub(':.*','', colnames(savercnt))
pmlist <- lapply(unique(allp), function(p){
  tmpmat <- savercnt[, allp == p]
  set.seed(1234567)
  colnames(tmpmat) <- sample(colnames(tmpmat))
  tmpmat
})
pmsavercnt <- do.call(cbind, pmlist)
pmsavercnt <- pmsavercnt[, pt[pt[,1] %in% colnames(savercnt), 1]] ## pseudotime for group1

pmlist <- lapply(unique(allp), function(p){
  tmpmat <- cnt[, allp == p]
  set.seed(1234567)
  colnames(tmpmat) <- sample(colnames(tmpmat))
  tmpmat
})
pmcnt <- do.call(cbind, pmlist)
pmcnt <- pmcnt[, pt[pt[,1] %in% colnames(cnt), 1]] ## pseudotime for group1

### permuted a group2 matrix for adding mean
allp2 <-  sub(':.*','', colnames(savercnt2))
pmlist2 <- lapply(unique(allp2), function(p){
  tmpmat <- savercnt2[, allp2 == p]
  set.seed(1234567)
  colnames(tmpmat) <- sample(colnames(tmpmat))
  tmpmat
})
pmsavercnt2 <- do.call(cbind, pmlist2)
pmsavercnt2 <- pmsavercnt2[, pt[pt[,1] %in% colnames(savercnt2), 1]] ## pseudotime for group2

allp2 <-  sub(':.*','', colnames(cnt2))
pmlist2 <- lapply(unique(allp2), function(p){
  tmpmat <- cnt2[, allp2 == p]
  set.seed(1234567)
  colnames(tmpmat) <- sample(colnames(tmpmat))
  tmpmat
})
pmcnt2 <- do.call(cbind, pmlist2)
pmcnt2 <- pmcnt2[, pt[pt[,1] %in% colnames(cnt2), 1]] ## pseudotime for group2

### demean a group2 matrix for adding trend
allp2 <-  sub(':.*','', colnames(savercnt2))
pmlist2 <- lapply(unique(allp2), function(p){
  tmpmat <- savercnt2[, allp2 == p]
  tmpmat - rowMeans(tmpmat)
})
savercnt2.demean <- do.call(cbind, pmlist2)
savercnt2.demean <- savercnt2.demean[, pt[pt[,1] %in% colnames(savercnt2), 1]] ## pseudotime for group2

allp2 <-  sub(':.*','', colnames(cnt2))
pmlist2 <- lapply(unique(allp2), function(p){
  tmpmat <- cnt2[, allp2 == p]
  tmpmat - rowMeans(tmpmat)
})
cnt2.demean <- do.call(cbind, pmlist2)
cnt2.demean <- cnt2.demean[, pt[pt[,1] %in% colnames(cnt2), 1]] ## pseudotime for group2

### add signal to  expression (both cnt and savercnt) 
dir.create('count/', showWarnings = F, recursive = T)
dir.create('saver/', showWarnings = F, recursive = T)
for (j in seq(1,4)) { # signal from 1 weakest to 4 strongest
  print(j)
  fromgene <- lapply(seq(1, max(clu)), function(i){
    print(i)
    clumat = saverlog[names(clu[clu==i]), ]  # log2 expression
    fstat.clu <- fstat[rownames(clumat), 1]
    s <- names(sort(fstat.clu))
    set.seed(1234567)
    fromgenetmp <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], floor(length(selgene)/10), replace = T)
    set.seed(1234567)
    sample(fromgenetmp)  ## randomize for three types of selgenes
  }) ## j = 1, smallest fstat, weakest signal
  fromgene = unlist(fromgene)
  if (length(fromgene) < length(selgene)){
    set.seed(1234567)
    fromgene = c(fromgene, sample(fromgene, length(selgene) - length(fromgene), replace = TRUE))
  }
  saveRDS(fromgene, paste0('fromgene/', j, '.rds'))
  
  ### SAVER
  m_with_halfmean <- pmsavercnt2[fromgene[((length(selgene1)+length(selgene2)+1):(length(selgene1)+length(selgene2)+length(selgene3)))],]
  m_with_halfmean <- m_with_halfmean/2
                            
  resexpr <-  rbind(savercnt[selgene1, ] + savercnt[fromgene[1:length(selgene1)], ], 
              savercnt[selgene2,] + pmsavercnt[fromgene[(length(selgene1)+1) : (length(selgene1)+length(selgene2))],],
              savercnt[selgene3,] + savercnt[fromgene[((length(selgene1)+length(selgene2)+1):(length(selgene1)+length(selgene2)+length(selgene3)))], ])
  mat <- rbind(resexpr[selgene, ], savercnt[setdiff(rownames(savercnt),selgene),])
  #### add three type of signals to group2 compared to group1: mean only, trend only, and mean and trend only.
  add1 <- rbind(savercnt2[selgene1, ] + pmsavercnt2[fromgene[1:length(selgene1)], ],
                savercnt2[selgene2, ],
                savercnt2[selgene3, ] + m_with_halfmean)
  savercnt2.new <- rbind(add1, savercnt2[setdiff(rownames(savercnt2), selgene), ]) ## others do not add
  #### bind group1 and group2
  mat <- cbind(mat, savercnt2.new[rownames(mat), ])
  saveRDS(mat, paste0('saver/', j,'.rds'))
  
  ### count
  m_with_halfmean <- pmcnt2[fromgene[((length(selgene1)+length(selgene2)+1):(length(selgene1)+length(selgene2)+length(selgene3)))],]
  m_with_halfmean <- m_with_halfmean/2
  
  resexpr <-  rbind(cnt[selgene1, ] + cnt[fromgene[1:length(selgene1)], ], 
              cnt[selgene2,] + pmcnt[fromgene[(length(selgene1)+1) : (length(selgene1)+length(selgene2))],],
              cnt[selgene3,] + cnt[fromgene[((length(selgene1)+length(selgene2)+1):(length(selgene1)+length(selgene2)+length(selgene3)))], ])
  mat <- rbind(resexpr[selgene, ], cnt[setdiff(rownames(cnt),selgene),])
  
  
  add1 <- rbind(cnt2[selgene1, ] + pmcnt2[fromgene[1:length(selgene1)], ],
                cnt2[selgene2, ],
                cnt2[selgene3, ] + m_with_halfmean)
  cnt2.new <- rbind(add1, cnt2[setdiff(rownames(cnt2), selgene), ])
  mat <- cbind(mat, cnt2.new[rownames(mat), ])
  saveRDS(mat, paste0('count/', j,'.rds'))
  
  ## double check
  # logmat <- log2(mat+1)
  # sum(logmat < 0 )
  # summary(sapply(sample(selgene,10),function(i) {  ## should be larger
  #   summary(lm(logmat[i,pt[pt[,1] %in% colnames(logmat),1]]~I(1:ncol(logmat))))$fstatistic[1]
  # }))
  # summary(sapply(sample(othgene,10),function(i) { ## should be much smaller
  #   summary(lm(logmat[i,pt[pt[,1] %in% colnames(logmat),1]]~I(1:ncol(logmat))))$fstatistic[1]
  # }))
}



