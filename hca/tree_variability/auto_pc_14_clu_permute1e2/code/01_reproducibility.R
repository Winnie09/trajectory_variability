rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
n.permute <- 100
suppressMessages(library(igraph))
setwd("/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate")

umap = readRDS('umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)

ctlevel <- data.frame(ct=c('HSC','MPP','LMPP','CMP','CLP','GMP','MEP',"Bcell","CD4Tcell","CD8Tcell",'NKcell','Mono','Ery'),level=c(1,2,3,3,4,4,4,5,5,5,5,5,5),immunepath=c(1,1,1,0,1,0,0,1,1,1,1,0,0),monopath=c(1,1,1,1,0,1,0,0,0,0,0,1,0),erypath=c(1,1,0,1,0,0,1,0,0,0,0,0,1),stringsAsFactors = F)
str(pca)

a = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/ct/sc.rds')
ct = data.frame(cell = names(a), celltype = a, stringsAsFactors = FALSE)


### determine numPC
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

set.seed(12345)
sdev <- apply(pca, 2, sd)
x <- 1:50
optpoint <- which.min(sapply(2:20, function(i) {
  x2 <- pmax(0, x - i)
  sum(lm(sdev ~ x + x2)$residuals^2)
}))
pcadim = optpoint + 1
pr <- pca[,1:pcadim]  # 7
# pr <- pca[,1:2]

## clustering
clu <- mykmeans(pr, number.cluster = 14)$cluster
pd = data.frame(x = pr[,1], y = pr[,2], clu = as.factor(clu[rownames(pr)]))

mypalette = colorRampPalette(brewer.pal(9,'Set1'))
ggplot(data = pd, aes(x = x, y = y, color = clu)) + 
  geom_scattermore()+
  scale_color_manual(values = mypalette(14))

## cell type composition in clusters
pd = cbind(pd, celltype = ct[match(rownames(pd), ct[,1]),2])
tab <- table(pd[,3:4])
tab <- tab/rowSums(tab)
pd <- melt(tab)
pd$clu <- factor(as.character(pd$clu), levels = seq(1,max(pd$clu)))

ggplot(data = pd) +
  geom_bar(aes(x = clu, y = value, fill = celltype), stat = 'identity', position = 'dodge') +
  theme_classic() +
  ylab('Celltype Proportion') +
  scale_fill_manual(values = mypalette(length(unique(pd$celltype))))

# tmp <- which.min(sapply(1:clun,function(scn) mean(ctlevel[match(ct[match(names(clu)[clu==scn],ct[,1]),3],ctlevel[,1]),2],na.rm=T)))
# 

###
mcl <- exprmclust(t(pr),cluster=clu,reduce=F)
# mcl <- exprmclust(t(pr), reduce = F)
plotmclust(mcl, cell_point_size = 0.1)
str(mcl)
## find origin
tmp <- pd[pd$celltype == 'HSC', ]
origin.cluster <- as.numeric(tmp[which.max(tmp[,3]), 1])


ord <- TSCANorder(mcl, startcluster = origin.cluster, listbranch = T,orderonly = T)
str(ord)
length(ord)
pt <- data.frame(cell = unname(unlist(ord)), time = c(1:length(ord[[1]]), 1:length(ord[[2]]), 1:length(ord[[3]]), 1:length(ord[[4]]), 1:length(ord[[5]]), 1:length(ord[[6]])), stringsAsFactors = FALSE)

## plot pseudotime
pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(pt[match(rownames(pca), pt[,1]),2]))
library(scattermore)
library(RColorBrewer)


ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
  geom_scattermore() +
  scale_color_gradient(low = 'yellow', high = 'blue')

# -----------
# permutation 
# -----------
jslist <- oclist <- list()
for (pmid in seq(1, n.permute)){
  ## boostrap cells
  print(pmid)
  set.seed(pmid)
  pr.pm <- pr[sample(rownames(pr), nrow(pr), replace = TRUE),]
  pr.pm <- pr.pm[!duplicated(rownames(pr.pm)),]
  
  ## cluster cells
  clu <- mykmeans(pr.pm, number.cluster = 14)$cluster ###
  mcl.pm <- exprmclust(t(pr.pm), cluster = clu, reduce = FALSE) ###
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
  library(scattermore)
  library(RColorBrewer)
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
str(jsm)
str(ocm)
par(mfrow = c(1,2))
hist(jsm)
hist(ocm)

js.cut <- 0.5
oc.cut <- 0.6

res <- sapply(seq(1,length(jslist)), function(i){
  js <- jslist[[i]]
  js.binary <- (js > js.cut) + 0

  while (length(which(rowSums(js.binary) > 1)) > 0 | length(which(colSums(js.binary) > 1)) > 0){
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
      
    dup.id <- which(colSums(js.binary) > 1)
    if (length(dup.id) == 1){
      addid <- which.max(js[, dup.id])
      js.binary[dup.id, ] <- 0
      js.binary[addid, dup.id] <- 1  
    } else if (length(dup.id) > 1) {
      for (dup.i in dup.id){
        addid <- which.max(js[, dup.id])
        js.binary[, dup.id] <- 0
        js.binary[addid, dup.id] <- 1  
      }
    }
  }
  js.melt <- melt(js.binary)
  js.melt <- js.melt[js.melt[,3]!=0,]
  as.character(js.melt[,2])
})
res <- unlist(res)  
js.perc <- table(res)/n.permute
saveRDS(js.perc, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/js_percentage.rds')

res <- sapply(seq(1,length(oclist)), function(i){
  oc <- oclist[[i]]
  oc.binary <- (oc > oc.cut) + 0
  while (length(which(rowSums(oc.binary) > 1)) > 0 | length(which(colSums(oc.binary) > 1)) > 0){
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
    
    dup.id <- which(colSums(oc.binary) > 1)
    if (length(dup.id) == 1){
      addid <- which.max(oc[, dup.id])
      oc.binary[, dup.id] <- 0
      oc.binary[addid, dup.id] <- 1  
    } else if (length(dup.id) > 1) {
      for (dup.i in dup.id){
        addid <- which.max(oc[, dup.i])
        oc.binary[, dup.i] <- 0
        oc.binary[addid, dup.i] <- 1  
      }
    }
  }
    
  oc.melt <- melt(oc.binary)
  oc.melt <- oc.melt[oc.melt[,3]!=0,]
  if (length(oc.melt[,2]) > 6) print(i)
  as.character(oc.melt[,2])
})
res <- unlist(res)  
oc.perc <- table(res)/n.permute
saveRDS(oc.perc, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/oc_percentage.rds')



