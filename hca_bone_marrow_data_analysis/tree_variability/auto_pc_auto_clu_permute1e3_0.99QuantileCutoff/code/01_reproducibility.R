rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
setwd("/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate")
umap = readRDS('umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)
ctlevel <- data.frame(ct=c('HSC','MPP','LMPP','CMP','CLP','GMP','MEP',"Bcell","CD4Tcell","CD8Tcell",'NKcell','Mono','Ery'),level=c(1,2,3,3,4,4,4,5,5,5,5,5,5),immunepath=c(1,1,1,0,1,0,0,1,1,1,1,0,0),monopath=c(1,1,1,1,0,1,0,0,0,0,0,1,0),erypath=c(1,1,0,1,0,0,1,0,0,0,0,0,1),stringsAsFactors = F)
str(pca)
a = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/ct/sc.rds')
ct = data.frame(cell = names(a), celltype = a, stringsAsFactors = FALSE)

mykmeans <- function(matrix, number.cluster = NA){
  ## cluster the rows
  set.seed(12345)
  library(parallel)
  if (is.na(number.cluster)){
    maxclunum <- 20
    rss <- mclapply(1:maxclunum,function(clunum) {
      tmp <- kmeans(matrix,clunum,iter.max = 1000)
      tmp$betweenss/tmp$totss
    },mc.cores=20)
    rss <- unlist(rss)
    x <- 1:maxclunum
    optclunum <- which.min(sapply(1:maxclunum, function(i) {
        x2 <- pmax(0, x - i)
        sum(lm(rss ~ x + x2)$residuals^2)  ## check this
    }))
    clu <- kmeans(matrix,optclunum)
  } else {
    clu <- kmeans(matrix, number.cluster)    
  }
    return(clu)
}

### determine numPC
set.seed(12345)
sdev <- apply(pca, 2, sd)
x <- 1:20
optpoint <- which.min(sapply(2:20, function(i) {
  x2 <- pmax(0, x - i)
  sum(lm(sdev[1:20] ~ x + x2)$residuals^2)
}))
pcadim = optpoint + 1
pr <- pca[,1:pcadim]  # 2

### mclust
# mcl <- exprmclust(t(pr),cluster=clu,reduce=F)
mcl <- exprmclust(t(pr), reduce = F)
plotmclust(mcl, cell_point_size = 0.1)
str(mcl)
## find origin
pd = data.frame(x = pr[,1], y = pr[,2], clu = as.factor(mcl$clusterid))
pd = cbind(pd, celltype = ct[match(rownames(pd), ct[,1]),2])
tab <- table(pd[,3:4])
tab <- tab/rowSums(tab)
pd <- melt(tab)
pd$clu <- factor(as.character(pd$clu), levels = seq(1,max(pd$clu)))

tmp <- pd[pd$celltype == 'HSC', ]
origin.cluster <- as.numeric(tmp[which.max(tmp[,3]), 1])

## construct pseudotime
ord <- TSCANorder(mcl, startcluster = origin.cluster, listbranch = T,orderonly = T)
str(ord)
length(ord)
pt <- data.frame(cell = unname(unlist(ord)), time = unlist(sapply(sapply(ord, length), function(i) seq(1, i))), stringsAsFactors = FALSE)

## plot pseudotime
pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(pt[match(rownames(pca), pt[,1]),2]))
library(scattermore)
library(RColorBrewer)
ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
  geom_scattermore() +
  scale_color_gradient(low = 'yellow', high = 'blue')

# -------------------------------------------------------
# null distribution of Jaccard index, overlap coefficient
# -------------------------------------------------------
js.null <- lapply(seq(1, length(ord)), function(i){
  b.ori <- ord[[i]]
  tmp <- sapply(seq(1, 1e3), function(j){
    set.seed(j)
    b.pm <- sample(rownames(pr), length(b.ori))
    length(intersect(b.pm, b.ori))/length(union(b.pm, b.ori))
  })
})

par(mfrow = c(1,3))
hist(js.null[[1]])
hist(js.null[[2]])
hist(js.null[[3]])

js.cut <- sapply(js.null, quantile, 0.99)

oc.null <- lapply(seq(1, length(ord)), function(i){
  b.ori <- ord[[i]]
  tmp <- sapply(seq(1, 1e3), function(j){
    set.seed(j)
    b.pm <- sample(rownames(pr), length(b.ori))
    length(intersect(b.pm, b.ori))/min(length(b.pm), length(b.ori))
  })
})
par(mfrow = c(1,3))
hist(oc.null[[1]])
hist(oc.null[[2]])
hist(oc.null[[3]])

oc.cut <- sapply(oc.null, quantile, 0.99)

# -----------
# permutation 
# -----------
jslist <- oclist <- list()
for (pmid in seq(1, 1e3)){
  ## boostrap cells
  print(pmid)
  set.seed(pmid)
  pr.pm <- pr[sample(rownames(pr), nrow(pr), replace = TRUE),]
  pr.pm <- pr.pm[!duplicated(rownames(pr.pm)),]
  
  # ## cluster cells
  mcl.pm <- exprmclust(t(pr.pm), reduce = FALSE) ###
  # plotmclust(mcl.pm, cell_point_size = 0.1)
  
  ## select origin cluster
  pt.pm.mean<- tapply(pt[match(names(mcl.pm[['clusterid']]), pt[,1]),2], list(mcl.pm[['clusterid']]), mean)
  start.cluster <- names(which.min(pt.pm.mean))
  
  ## construct pseudotime
  ord.pm <- TSCANorder(mcl.pm, startcluster = start.cluster, listbranch = T,orderonly = T)
  # str(ord.pm)
  
  ## plot pseudotime
  
  pt.pm <- data.frame(cell = unname(unlist(ord.pm)), time = unlist(sapply(sapply(ord.pm, length), function(i) seq(1, i))), stringsAsFactors = FALSE)
  pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(pt.pm[match(rownames(pca), pt.pm[,1]),2]))
  # ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
  #   geom_scattermore()

  ## compare two MST
  js <- sapply(seq(1, length(ord)), function(i){
          sapply(seq(1, length(ord.pm)), function(j){
            b.ori <- ord[[i]]
            b.pm <- ord.pm[[j]]
            js <- length(intersect(b.pm, b.ori))/length(union(b.pm, b.ori))
          })
      })
  
  oc <- sapply(seq(1, length(ord)), function(i){
           sapply(seq(1, length(ord.pm)), function(j){
              b.ori <- ord[[i]]
              b.pm <- ord.pm[[j]]
              oc <- length(intersect(b.pm, b.ori))/min(length(b.pm), length(b.ori))
           }) 
        })
  colnames(oc) <- colnames(js) <- paste0('original', seq(1, length(ord)))
  jslist[[pmid]] <- js
  oclist[[pmid]] <- oc
}
saveRDS(jslist, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/pm_js.rds')   
saveRDS(oclist, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/pm_oc.rds')   

jsm <- do.call(rbind, jslist)
ocm <- do.call(rbind, oclist)
par(mfrow = c(1,2))
hist(jsm)
hist(ocm)

res <- sapply(seq(1,length(jslist)), function(i){
  js <- jslist[[i]]
  js.binary <- sapply(seq(1,ncol(js)), function(c){
    (js[,c] > js.cut[c]) + 0
  })
  dup.id <- which(rowSums(js.binary) > 1)
  if (length(dup.id) == 1){
    addid <- which.max(js[dup.id, ])
    js.binary[dup.id, ] <- 0
    js.binary[dup.id, addid] <- 1  
  } else if (length(dup.id) > 1) {
    for (dup.i in dup.id){
      addid <- which.max(js[dup.i, ])
      js.binary[dup.i, ] <- 0
      js.binary[dup.i, addid] <- 1  
    }
  }
  js.melt <- melt(js.binary)
  js.melt <- js.melt[js.melt[,3]!=0,]
  as.character(js.melt[,2])
})
res <- unlist(res)  
js.perc <- table(res)/1e3
saveRDS(js.perc, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/js_percentage.rds')

res <- sapply(seq(1,length(oclist)), function(i){
  oc <- oclist[[i]]
  oc.binary <- sapply(seq(1,ncol(oc)), function(c){
    (oc[,c] > oc.cut[c]) + 0
  })
  dup.id <- which(rowSums(oc.binary) > 1)
  if (length(dup.id) == 1){
    addid <- which.max(oc[dup.id, ])
    oc.binary[dup.id, ] <- 0
    oc.binary[dup.id, addid] <- 1  
  } else if (length(dup.id) > 1) {
    for (dup.i in dup.id){
      addid <- which.max(oc[dup.i, ])
      oc.binary[dup.i, ] <- 0
      oc.binary[dup.i, addid] <- 1  
    }
  }
  oc.melt <- melt(oc.binary)
  oc.melt <- oc.melt[oc.melt[,3]!=0,]
  as.character(oc.melt[,2])
})
res <- unlist(res)  
oc.perc <- table(res)/1e3
sort((js.perc + oc.perc)/2)
saveRDS(oc.perc, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/oc_percentage.rds')


