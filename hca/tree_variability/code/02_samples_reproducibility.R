rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
n.permute <- 1e3
max.clunum <- 50
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
x <- 1:max.clunum
optpoint <- which.min(sapply(2:max.clunum, function(i) {
  x2 <- pmax(0, x - i)
  sum(lm(sdev[1:max.clunum] ~ x + x2)$residuals^2)
}))
pcadim = optpoint + 1
pr <- pca[,1:pcadim]  # 7

## clustering
clu <- mykmeans(pr, number.cluster = 14)$cluster
pd = data.frame(x = pr[,1], y = pr[,2], clu = as.factor(clu[rownames(pr)]))

# mypalette = colorRampPalette(brewer.pal(9,'Set1'))
# ggplot(data = pd, aes(x = x, y = y, color = clu)) + 
#   geom_scattermore()+
#   scale_color_manual(values = mypalette(14))+
#   theme_classic() + xlab('UMAP1') + ylab('UMAP2')

## cell type composition in clusters
pd = cbind(pd, celltype = ct[match(rownames(pd), ct[,1]),2])
tab <- table(pd[,3:4])
tab <- tab/rowSums(tab)
pd <- melt(tab)
pd$clu <- factor(as.character(pd$clu), levels = seq(1,max(pd$clu)))

# ggplot(data = pd) +
#   geom_bar(aes(x = clu, y = value, fill = celltype), stat = 'identity', position = 'dodge') +
#   theme_classic() +
#   ylab('Celltype Proportion') +
#   scale_fill_manual(values = mypalette(length(unique(pd$celltype))))

### mclust
mcl <- exprmclust(t(pr),cluster=clu,reduce=F)
# mcl <- exprmclust(t(pr), reduce = F)
# plotmclust(mcl, cell_point_size = 0.1)
# str(mcl)

# --------------------
# construct pseudotime 
# --------------------
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
pt <- unlist(sapply(sapply(ord, length), function(i) seq(1, i)))
names(pt) <- unname(unlist(ord))

# ## plot pseudotime
# pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(pt[rownames(pca)]))
# library(scattermore)
# library(RColorBrewer)
# ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
#   geom_scattermore() +
#   scale_color_gradient(low = 'yellow', high = 'blue')

# ------------------------------------------------------------
# get candidate branches to test reproducibility, 20200726 >>
# ------------------------------------------------------------
findbranch <- function(mst, order, origin){
  deg <- degree(mst)
  vertex <- names(deg[which(deg > 2 | deg == 1)])
  if (!origin %in% vertex) vertex <- c(origin, vertex)
  eg <- expand.grid(1:length(vertex), 1:length(vertex))
  eg <- eg[eg[,1]<eg[,2],]
  eg = data.frame(vertex[eg[,1]], vertex[eg[,2]], stringsAsFactors = FALSE)
  library(igraph)
  tmpbranch <- lapply(seq(1,nrow(eg)), function(i){
    sp <- shortest_paths(mst, from = eg[i,1], to = eg[i,2])$vpath[[1]]
    if (sum(vertex %in% sp) == 2) as.vector(sp)
  })
  tmpbranch <- tmpbranch[sapply(tmpbranch, length) >0]  
 
  allbranch <- gsub('backbone ', '', gsub('branch: ', '', names(order)))
  allbranch <- sapply(allbranch, function(i) strsplit(i, ',')[[1]])
  allbranch <- paste0(names(allbranch), collapse = ' ')
  newbranch <-sapply(tmpbranch, function(i) {
      tmp <- paste0(i, collapse = ',')
      if (!grepl(tmp, allbranch)){
        rev(i)
      } else {
        i
      }
  })
  return(newbranch)
}
newbranch <- findbranch(mst = mcl$MSTtree, order = ord, origin = origin.cluster)  


# -------------------------------------------------------
# null distribution of Jaccard index, overlap coefficient
# -------------------------------------------------------
## add here --------------->>>>>>
## for samples
## add here ---------------<<<<<<<
js.null <- lapply(seq(1, length(newbranch)), function(i){
  b.ori <- unlist(sapply(newbranch[[i]], function(c) names(mcl$clusterid[mcl$clusterid == c])))
  b.ori.alls <- gsub(':.*', '', b.ori)
  alls <- gsub(':.*', '', rownames(pr))
  tmp <- mclapply(seq(1, 1e3), function(j){
    set.seed(j)
    b.pm <- sample(rownames(pr), length(b.ori))
    b.pm.alls <- gsub(':.*', '', b.pm)
    tmpp <- sapply(unique(alls), function(s){
      b.pm.s <- b.pm[b.pm.alls == s]
      b.ori.s <- b.ori[b.ori.alls == s]
      length(intersect(b.pm.s, b.ori.s))/length(union(b.pm.s, b.ori.s))
    })  
  },mc.cores = detectCores()-2)
  tmp <- do.call(rbind,tmp)
})
js.cut <- sapply(js.null, function(i) apply(i, 2, quantile, 0.99))
# ------------------

oc.null <- lapply(seq(1, length(newbranch)), function(i){
  b.ori <- unlist(sapply(newbranch[[i]], function(c) names(mcl$clusterid[mcl$clusterid == c])))
  b.ori.alls <- gsub(':.*', '', b.ori)
  tmp <- mclapply(seq(1, 1e3), function(j){
    set.seed(j)
    b.pm <- sample(rownames(pr), length(b.ori))
    b.pm.alls <- gsub(':.*', '', b.pm)
    tmpp <- sapply(unique(alls), function(s){
      b.pm.s <- b.pm[b.pm.alls == s]
      b.ori.s <- b.ori[b.ori.alls == s]
      length(intersect(b.pm.s, b.ori.s))/min(length(b.pm.s), length(b.ori.s))
    })  
  },mc.cores = detectCores()-2)
  tmp <- do.call(rbind,tmp)
})
oc.cut <- sapply(oc.null, function(i) apply(i, 2, quantile, 0.99))

# -----------
# permutation 
# -----------
corrlist.alls <- jslist.alls <- oclist.alls <- list()
n.permute = 100
for (pmid in seq(1, n.permute)){
  ## boostrap cells
  print(pmid)
  set.seed(pmid)
  pr.pm <- pr[sample(rownames(pr), nrow(pr), replace = TRUE),]
  pr.pm <- pr.pm[!duplicated(rownames(pr.pm)),]
  
  # ## cluster cells
  clu <- mykmeans(pr.pm, number.cluster = 14)$cluster ###

  pd = data.frame(x = pr[names(clu),1], y = pr[names(clu),2], clu = as.factor(clu))
  pd.text.x = tapply(pd[,1], list(pd$clu), mean)
  pd.text.y = tapply(pd[,2], list(pd$clu), mean)
  pd.text = data.frame(x = pd.text.x, y = pd.text.y, clu = names(pd.text.x))
  pd.text[14,1:2] =  c(pd.text[14,1] + 2, pd.text[14,2] + 1)

  # ggplot() + 
  #   geom_scattermore(data = pd, aes(x = x, y = y, color = clu))+
  #   scale_color_manual(values = mypalette(14))+
  #   theme_classic() + xlab('UMAP1') + ylab('UMAP2') +
  #   geom_text(data = pd.text, aes(x = x, y = y, label = clu))
  

  ## cell type composition in clusters
  # pd = cbind(pd, celltype = ct[match(rownames(pd), ct[,1]),2])
  # tab <- table(pd[,3:4])
  # tab <- tab/rowSums(tab)
  # pd <- melt(tab)
  # pd$clu <- factor(as.character(pd$clu), levels = seq(1,max(pd$clu)))
  # ggplot(data = pd) +
  #   geom_bar(aes(x = clu, y = value, fill = celltype), stat = 'identity', position = 'dodge') +
  #   theme_classic() +
  #   ylab('Celltype Proportion') +
  #   scale_fill_manual(values = mypalette(length(unique(pd$celltype))))

  # build pseudotime
  mcl.pm <- exprmclust(t(pr.pm), cluster = clu, reduce = FALSE) ###
  # plotmclust(mcl.pm, cell_point_size = 0.1)
  
  ## select origin cluster
  pt.pm.mean<- tapply(pt[names(mcl.pm[['clusterid']])], list(mcl.pm[['clusterid']]), mean)
  start.cluster <- names(which.min(pt.pm.mean))
  
  ## construct pseudotime
  ord.pm <- TSCANorder(mcl.pm, startcluster = start.cluster, listbranch = T,orderonly = T)
  # str(ord.pm)
  
  ## plot pseudotime
  pt.pm <- unlist(sapply(sapply(ord.pm, length), function(i) seq(1, i)))
  names(pt.pm) <- unname(unlist(ord.pm))
  pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(pt.pm[rownames(pca)]))
  # ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
  #   geom_scattermore() + theme_classic()
  
  # get candidate branches
  newbranch.pm <- findbranch(mst = mcl.pm$MSTtree, order = ord.pm, origin = start.cluster)
  
  ## compare two MST
  js <- sapply(seq(1, length(newbranch)), function(i){
          print('i')
          print(i)
          id <- which(sapply(paste0(names(ord),','), function(k) grepl(paste0(paste0(newbranch[[i]], collapse = ','),','), k)))[1]
          cells <- ord[[id]]
          b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k) names(mcl$clusterid)[mcl$clusterid == k])), cells)
          b.ori.alls <- gsub(':.*', '', b.ori)
          tmp <- mclapply(seq(1, length(newbranch.pm)), function(j){
            print(j)
            id <- which(sapply(paste0(names(ord.pm),','), function(k) grepl(paste0(paste0(newbranch.pm[[j]], collapse = ','),','), k)))[1]
            cells <- ord.pm[[id]]
            b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k) names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
            b.pm.alls <- gsub(':.*', '', b.pm)
            # js <- length(intersect(b.pm, b.ori))/length(union(b.pm, b.ori))
            tmpp <- sapply(unique(alls), function(s){
              b.pm.s <- b.pm[b.pm.alls == s]
              b.ori.s <- b.ori[b.ori.alls == s]
              length(intersect(b.pm.s, b.ori.s))/length(union(b.pm.s, b.ori.s))
            })  
          },mc.cores = detectCores()-2)
          tmp <- do.call(rbind,tmp)
          rownames(tmp) <- paste0('branch.pm', seq(1, length(newbranch.pm)))
          tmp
        }, simplify = FALSE)
  names(js) <- paste0('branch', seq(1, length(newbranch)))
  ###### =====================================
  oc <- sapply(seq(1, length(newbranch)), function(i){
            id <- which(sapply(paste0(names(ord),','), function(k) grepl(paste0(paste0(newbranch[[i]], collapse = ','),','), k)))[1]
            cells <- ord[[id]]
            b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k) names(mcl$clusterid)[mcl$clusterid == k])), cells)
            b.ori.alls <- gsub(':.*', '', b.ori)
            tmp <- mclapply(seq(1, length(newbranch.pm)), function(j){
                id <- which(sapply(paste0(names(ord.pm),','), function(k) grepl(paste0(paste0(newbranch.pm[[j]], collapse = ','),','), k)))[1]
                cells <- ord.pm[[id]]
                b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k) names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
                b.pm.alls <- gsub(':.*', '', b.pm)
                # oc <- length(intersect(b.pm, b.ori))/min(length(b.pm), length(b.ori))
                tmpp <- sapply(unique(alls), function(s){
                  b.pm.s <- b.pm[b.pm.alls == s]
                  b.ori.s <- b.ori[b.ori.alls == s]
                  length(intersect(b.pm.s, b.ori.s))/min(length(b.pm.s), length(b.ori.s))
                })  
          },mc.cores = detectCores()-2)
          tmp <- do.call(rbind,tmp)
          rownames(tmp) <- paste0('branch.pm', seq(1, length(newbranch.pm)))
          tmp
        }, simplify = FALSE)
  names(oc) <- paste0('branch', seq(1, length(newbranch)))           
        
  
  corr <- sapply(seq(1, length(newbranch)), function(i){
              id <- which(sapply(paste0(names(ord),','), function(k) grepl(paste0(paste0(newbranch[[i]], collapse = ','),','), k)))[1]
              cells <- ord[[id]]
              b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k) names(mcl$clusterid)[mcl$clusterid == k])), cells)
              b.ori.alls <- gsub(':.*', '', b.ori)
              tmp <- mclapply(seq(1, length(newbranch.pm)), function(j){
                  id <- which(sapply(paste0(names(ord.pm),','), function(k) grepl(paste0(paste0(newbranch.pm[[j]], collapse = ','),','), k)))[1]
                  cells <- ord.pm[[id]]
                  b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k) names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
                  b.pm.alls <- gsub(':.*', '', b.pm)
                  # ov = intersect(b.ori, b.pm)
                  # cor(pt[ov], pt.pm[ov])
                  tmpp <- sapply(unique(alls), function(s){
                    b.pm.s <- b.pm[b.pm.alls == s]
                    b.ori.s <- b.ori[b.ori.alls == s] 
                    ov = intersect(b.ori.s, b.pm.s)
                    cor(pt[ov], pt.pm[ov])
                  })
              }, mc.cores = detectCores()-2) 
              tmp <- do.call(rbind, tmp)
              rownames(tmp) <- paste0('branch.pm', seq(1, length(newbranch.pm)))
              tmp[is.na(tmp)] <- 0
              tmp
          }, simplify = FALSE)
  # corr[is.na(corr)] <- 0
  names(corr) <- paste0('branch', seq(1, length(newbranch)))           
  # colnames(corr) <- colnames(oc) <- colnames(js) <- paste0('original', seq(1, length(newbranch)))
  jslist.alls[[pmid]] <- js
  oclist.alls[[pmid]] <- oc
  corrlist.alls[[pmid]] <- corr
}
saveRDS(jslist.alls, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/samples/pm_js_alls.rds')   
saveRDS(oclist.alls, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/samples/pm_oc_alls.rds')   
saveRDS(corrlist.alls, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/samples/result/pm_oc_alls.rds')   

# jsm <- do.call(rbind, jslist)
# ocm <- do.call(rbind, oclist)
# par(mfrow = c(1,2))
# hist(jsm)
# hist(ocm)
s = unique(alls)[1]
df.alls <- lapply(unique(alls), function(s){
  jslist = sapply(jslist.alls, function(i){
    sapply(i, function(ii) ii[,s])
  }, simplify = FALSE)
  oclist = sapply(oclist.alls, function(i){
    sapply(i, function(ii) ii[,s])
  }, simplify = FALSE)
  corrlist = sapply(corrlist.alls, function(i){
    sapply(i, function(ii) ii[,s])
  }, simplify = FALSE)
  
  res <- corr.score <- list()
  for (i in seq(1, length(jslist))){
    print(i)
    js <- jslist[[i]]
    js.binary <- sapply(seq(1,ncol(js)), function(c){
      (js[,c] > js.cut[c]) + 0
    })
    while (length(which(rowSums(js.binary) > 1)) > 0 | length(which(colSums(js.binary) > 1)) > 0){
      dup.id <- which(rowSums(js.binary) > 1)
      if (length(dup.id) == 1){
        addid <- which.max(js[dup.id, ])
        js.binary[dup.id, ] <- 0
        js.binary[dup.id, addid] <- 1  
      } else if (length(dup.id) > 1) {
        for (dup.i in dup.id){
          print(dup.i)
          addid <- which.max(js[dup.i, ])
          js.binary[dup.i, ] <- 0
          js.binary[dup.i, addid] <- 1  
        }
      }
        
      dup.id <- which(colSums(js.binary) > 1)
      if (length(dup.id) == 1){
        addid <- which.max(js[, dup.id])
        js.binary[, dup.id] <- 0
        js.binary[addid, dup.id] <- 1  
      } else if (length(dup.id) > 1) {
        for (dup.i in dup.id){
          addid <- which.max(js[, dup.i])
          js.binary[, dup.i] <- 0
          js.binary[addid, dup.i] <- 1  
        }
      }
    }
    
    corr.score[[i]] <- corrlist[[i]] * js.binary
    js.melt <- melt(js.binary)
    js.melt <- js.melt[js.melt[,3]!=0,]
    res[[i]] <- as.character(js.melt[,2])
  }
  res <- unlist(res)  
  js.perc <- table(res)/n.permute
  names(js.perc) <- newbranch
  # saveRDS(js.perc, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/samples/js_percentage.rds')
  
  corr.score.m <- do.call(rbind, corr.score)
  corr.score.v <- colSums(corr.score.m)/n.permute
  names(corr.score.v) <- newbranch
  # saveRDS(corr.score.v, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/samples/corr_score.rds')
  
  res <- sapply(seq(1,length(oclist)), function(i){
    print(i)
    oc <- oclist[[i]]
    oc.binary <- sapply(seq(1,ncol(oc)), function(c){
      (oc[,c] > oc.cut[c]) + 0
    })
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
    as.character(oc.melt[,2])
  })
  res <- unlist(res)  
  oc.perc <- table(res)/n.permute
  names(oc.perc) <- newbranch
  sort((js.perc + oc.perc)/2)
  
  df <- data.frame(js.perc = js.perc, oc.perc = oc.perc, corr.score.v = corr.score.v)
  df <- df[, c(2,4,5)]

# saveRDS(oc.perc, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/result/samples/oc_percentage.rds')
  
})
names(df.alls) <- unique(alls)
df.alls[order(names(df.alls))]



