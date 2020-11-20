## 03_preparedata.R
geneProp <- 0.2
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/poolSampleSignal/')
### load saver, and count matrix
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
saverlog <- readRDS('./null/hsc_mep_ery_saver.rds')
cnt <- readRDS('./null/hsc_mep_ery_count.rds')

### prepare count, imputed, selected genes
savercnt <- 2^saverlog - 1
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
rownames(savercnt) <- sapply(rownames(savercnt), function(i) sub('_','-',i))
cnt <- cnt[, colnames(savercnt)]
allp = sub(':.*','', colnames(savercnt))
sample <- sub(':.*','',colnames(savercnt))
names(sample) <- colnames(savercnt)
set.seed(12345)
selgene <- sample(row.names(savercnt), round(geneProp * nrow(savercnt)))
othgene <- setdiff(rownames(savercnt), selgene)
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

# ========  
### order cells by pseudotime, and get the optimal clusters using saver imputed 80% genes
pt <- readRDS('./null/pseudotime.rds')
savercnt <- savercnt[, pt[,1]]
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

# consider integrated all samples on pt, fit a curve along pt, use this curve's sd as the signal strength
    # predict the expr on pt
    # calculate the sd for this fitting curve\
    # order othgene (non-goldstandard) by sd, low to high
    # add signal of the averaged curve to each sample on goldstandard genes
    #In this case, the signals added to gs genes will have same pattern for all samples. The background noise should be kept the same.
fitmat <- get_spline_fit(trainData = savercnt[othgene,pt[,1]], trainX = pt[,2], fit.min = min(pt[,2]), fit.max = max(pt[,2]), fit.num.points = nrow(pt))
colnames(fitmat) <- colnames(savercnt)    
sds <- apply(fitmat, 1, sd)
# g = "DPM1:ENSG00000000419" 
# smoothScatter(savercnt[g,pt[,1]]~pt[,2],pch=20, col='grey',main=g,xlab='pseudotime',ylab='SAVER-imputed count')
# points(fitmat[g,pt[,1]]~pt[,2], pch = 20, col='red')


fitcntmat <- get_spline_fit(trainData = cnt[othgene,pt[,1]], trainX = pt[,2], fit.min = min(pt[,2]), fit.max = max(pt[,2]), fit.num.points = nrow(pt))
colnames(fitcntmat) <- colnames(savercnt)   
# g = "DPM1:ENSG00000000419" 
# smoothScatter(cnt[g,pt[,1]]~pt[,2],pch=20, col='grey',main=g,xlab='pseudotime',ylab='count')
# points(fitcntmat[g,pt[,1]]~pt[,2], pch = 20, col='red')

for( i in seq(1,10)){
  print(i)
  for (j in seq(1,4)) {
    # i = 9
    # fitclumat <- fitmat[names(clu[clu==i]),]

    ## add signals of othgene to selgene !!!
    
    s <- names(sort(sds[names(clu[clu == i])], decreasing = TRUE))
    
    # j = 1:4
    set.seed(12345)
    addgene <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], length(selgene), replace = T)
    addsavercnt <- fitmat[addgene, pt[,1]]
    ressavercnt <-  pmsavercnt[selgene, pt[,1]] + addsavercnt
    
    mat <- rbind(ressavercnt[selgene, pt[,1]], pmsavercnt[othgene, pt[,1]])
    saveRDS(mat, paste0('saver/clusterType', i, '_', j,'.rds'))
    saveRDS(list(selgene = selgene, addgene = rownames(addsavercnt)), paste0('saver/clusterType', i, '_', j,'_selgene_addgene.rds'))
    # this way to generate count data might be inaccurately
    addcnt <- fitcntmat[addgene, pt[,1]]
    rescnt <-  pmcnt[selgene, pt[,1]] + addcnt
    mat <- rbind(rescnt[selgene, pt[,1]], pmcnt[othgene, pt[,1]])
    tmp <- matrix(pmax(0, mat), nrow = nrow(mat))
    dimnames(tmp) <- dimnames(mat)
    saveRDS(tmp, paste0('count/clusterType', i, '_', j,'.rds'))
  }
}
    
# -------------------------------
# # inserted: plot example genes
# g <- "RPL39:ENSG00000198918" 
# g <- rownames(cnt)[grepl('^DHX37:', rownames(cnt))]
# library(ggplot2)
# library(reshape2)
# pd1 = data.frame(cnt = pmsavercnt[g, pt[,1]], 
#                 pseudotime = pt[,2], type = 'pmsavercnt', sample = sub(':.*', '', pt[,1]))
# pd2 = data.frame(cnt = addsavercnt[which(selgene == g), pt[,1]],
#                 pseudotime = pt[,2], type = 'addsavercnt', sample = sub(':.*', '', pt[,1]))
# pd3 = data.frame(cnt = ressavercnt[g, pt[,1]],
#                 pseudotime = pt[,2], type = 'ressavercnt', sample = sub(':.*', '', pt[,1]))
# ggplot(data = rbind(pd1, pd2, pd3)) + geom_point(aes(y = cnt, x = pseudotime, color = sample), size = 0.1) +
#   theme_classic() +
#   facet_wrap(~type) +
#   xlab('pseudotime') + ylab('savercount') +
#   theme(axis.text.x = element_blank())
# 
# pd = rbind(pd1, pd2, pd3)
# pd = pd[order(pd$sample, pd$pseudotime), ]
# pd$order = rep(1:nrow(pd1), each = 3)
# ggplot(data = pd) + geom_point(aes(y = cnt, x = order, color = sample, group = sample), size = 0.1) +
#   theme_classic() +
#   facet_wrap(~type) +
#   xlab('pseudotime(stratify sample)') + ylab('savercount') +
#   theme(axis.text.x = element_blank())


# ------------------------------
# # plot
# par(mfrow=c(1,3))
# plot(pt[,2],pmsavercnt['HBB:ENSG00000244734', pt[,1]], cex=0.1) 
# plot(pt[,2],savercnt['HBB:ENSG00000244734', pt[,1]], cex=0.1) 
# plot(pt[,2],ressavercnt['HBB:ENSG00000244734', pt[,1]], col='red', cex=0.1)
# 
# par(mfrow=c(1,3))
# plot(pt[,2],pmsavercnt[selgene[2], pt[,1]], cex=0.1) 
# plot(pt[,2],savercnt[addgene[2], pt[,1]], cex=0.1) 
# plot(pt[,2],ressavercnt[selgene[2], pt[,1]], col='red', cex=0.1)
# 
# 
# cv1 <- apply(ressavercnt[selgene, ], 1, sd)/rowMeans(ressavercnt[selgene, ])
# cv2 <- apply(pmsavercnt[othgene, ], 1, sd)/rowMeans(pmsavercnt[othgene, ])
# ggplot(data=rbind(data.frame(cv = cv1, g = 'selgene'), data.frame(cv = cv2, g = 'othgene')), aes(x=g, y = cv, col=g)) + geom_boxplot(outlier.colour = 'NA')
# 
# c1 <- sapply(selgene, function(g) cor(ressavercnt[g, ], seq(1, ncol(ressavercnt))))
# c2 <- sapply(othgene, function(g) cor(pmsavercnt[g, ], seq(1, ncol(pmsavercnt))))
# ggplot(data=rbind(data.frame(c = c1, g = 'selgene'), data.frame(c = c2, g = 'othgene')), aes(x=g, y = c, col=g)) + geom_boxplot(outlier.colour = 'NA')
# 
# 
# mat <- rbind(ressavercnt[selgene, pt[,1]], pmsavercnt[othgene, pt[,1]])
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

