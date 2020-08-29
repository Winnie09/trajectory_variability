rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
n.permute <- 3
max.clunum <- 50
setwd("/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate")

# --------------------------------------------------------------
# input: seurat integrated object including:
#  umap, pca
# celltype: a dataframe, col 1 is cell name, col 2 is cell type
# origin: the origin cell type
# --------------------------------------------------------------
# read in data
umap = readRDS('umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)
# ctlevel <- data.frame(ct=c('HSC','MPP','LMPP','CMP','CLP','GMP','MEP',"Bcell","CD4Tcell","CD8Tcell",'NKcell','Mono','Ery'),level=c(1,2,3,3,4,4,4,5,5,5,5,5,5),immunepath=c(1,1,1,0,1,0,0,1,1,1,1,0,0),monopath=c(1,1,1,1,0,1,0,0,0,0,0,1,0),erypath=c(1,1,0,1,0,0,1,0,0,0,0,0,1),stringsAsFactors = F)
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
get_binary <- function(matrix, matrix.cut){
  ## match boostrap and origin branches.
  ## matrix: #boostrap.branch * #origin.branch, values are js or oc
  ## matrix.cut: js or oc null distribution cutoff 
  matrix.binary <- sapply(seq(1,ncol(matrix)), function(c){
    (matrix[,c] > matrix.cut[c]) + 0
  })
  while (length(which(rowSums(matrix.binary) > 1)) > 0 | length(which(colSums(matrix.binary) > 1)) > 0){
    dup.id <- which(rowSums(matrix.binary) > 1)
    if (length(dup.id) == 1){
      addid <- which.max(matrix[dup.id, ])
      matrix.binary[dup.id, ] <- 0
      matrix.binary[dup.id, addid] <- 1  
    } else if (length(dup.id) > 1) {
      for (dup.i in dup.id){
        addid <- which.max(matrix[dup.i, ])
        matrix.binary[dup.i, ] <- 0
        matrix.binary[dup.i, addid] <- 1  
      }
    }
      
    dup.id <- which(colSums(matrix.binary) > 1)
    if (length(dup.id) == 1){
      addid <- which.max(matrix[, dup.id])
      matrix.binary[, dup.id] <- 0
      matrix.binary[addid, dup.id] <- 1  
    } else if (length(dup.id) > 1) {
      for (dup.i in dup.id){
        addid <- which.max(matrix[, dup.i])
        matrix.binary[, dup.i] <- 0
        matrix.binary[addid, dup.i] <- 1  
      }
    }
  }
  return(matrix.binary)
}

### determine numPC
infer_tree_structure <- function(pca, ct, origin.celltype){
  alls <- sub(':.*', '', ct$cell)
  names(alls) <- ct$cell
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
  # pd = data.frame(x = pr[,1], y = pr[,2], clu = as.factor(clu[rownames(pr)]))
  # mypalette = colorRampPalette(brewer.pal(9,'Set1'))
  # ggplot(data = pd, aes(x = x, y = y, color = clu)) + 
  #   geom_scattermore()+
  #   scale_color_manual(values = mypalette(14))+
  #   theme_classic() + xlab('UMAP1') + ylab('UMAP2')
  
  # ## cell type composition in clusters
  # pd = cbind(pd, celltype = ct[match(rownames(pd), ct[,1]),2])
  # tab <- table(pd[,3:4])
  # tab <- tab/rowSums(tab)
  # pd <- melt(tab)
  # pd$clu <- factor(as.character(pd$clu), levels = seq(1,max(pd$clu)))
  # 
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
  tmp <- pd[pd$celltype == origin.celltype, ]
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
  
  newbranch <- findbranch(mst = mcl$MSTtree, order = ord, origin = origin.cluster)  
  
  # -----------------------------------------------------
  # Evaluate robustness of tree branches using resampling
  # -----------------------------------------------------
  
  # null distribution of Jaccard index, overlap coefficient
  
  js.null <- lapply(seq(1, length(newbranch)), function(i) {
    b.ori <- unlist(sapply(newbranch[[i]], function(c) names(mcl$clusterid[mcl$clusterid == c])))
    tmp <- sapply(seq(1, 1e3), function(j){
      set.seed(j)
      b.pm <- sample(rownames(pr), length(b.ori))
      length(intersect(b.pm, b.ori))/length(union(b.pm, b.ori))
    })
  })
  
  # par(mfrow = c(2,ceiling(length(js.null)/2)))
  # for (i in js.null) hist(i, xlab = 'js', main = '', breaks = 50)
  
  js.cut <- sapply(js.null, quantile, 0.99)
  
  oc.null <- lapply(seq(1, length(newbranch)), function(i){
    b.ori <- unlist(sapply(newbranch[[i]], function(c) names(mcl$clusterid[mcl$clusterid == c])))
    tmp <- sapply(seq(1, 1e3), function(j){
      set.seed(j)
      b.pm <- sample(rownames(pr), length(b.ori))
      length(intersect(b.pm, b.ori))/min(length(b.pm), length(b.ori))
    })
  })
  # par(mfrow = c(2,ceiling(length(oc.null)/2)))
  # for (i in oc.null) hist(i, xlab = 'oc', main = '', breaks = 50)
  oc.cut <- sapply(oc.null, quantile, 0.99)
  
  mcl$pseudotime <- pt
  mcl$branch <- newbranch
  mcl$js.cut <- js.cut
  mcl$oc.cut <- oc.cut
  mcl$pca <- pr
  mcl$order <- ord
  mcl$allsample <- alls
  return(mcl)
}
# permutation 

a = infer_tree_structure(pca = pca, ct = ct, origin.celltype = 'HSC')

evaluate_uncertainty <- function(inferobj, n.permute){
  pr <- inferobj$pca
  newbranch <- inferobj$branch
  js.cut <- inferobj$js.cut
  oc.cut <- inferobj$oc.cut 
  pt <- inferobj$pseudotime
  ord <- inferobj$order
  alls <- inferobj$allsample
  ctcomplist <- reproduce.js <- reproduce.oc <- corr.score <- list()
  for (pmid in seq(1, n.permute)){
    ## boostrap cells
    print(pmid)
    set.seed(pmid)
    pr.pm <- pr[sample(rownames(pr), nrow(pr), replace = TRUE),]
    pr.pm <- pr.pm[!duplicated(rownames(pr.pm)),]
    
    ## cluster cells
    clu <- mykmeans(pr.pm, number.cluster = 14)$cluster ###
  
    # --- check if these codes are necessary <<<<<<<<<<<<<<<<
    # pd = data.frame(x = pr[names(clu),1], y = pr[names(clu),2], clu = as.factor(clu))
    # pd.text.x = tapply(pd[,1], list(pd$clu), mean)
    # pd.text.y = tapply(pd[,2], list(pd$clu), mean)
    # pd.text = data.frame(x = pd.text.x, y = pd.text.y, clu = names(pd.text.x))
    # pd.text[14,1:2] =  c(pd.text[14,1] + 2, pd.text[14,2] + 1)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
  
    ## build pseudotime
    mcl.pm <- exprmclust(t(pr.pm), cluster = clu, reduce = FALSE) ###
    # plotmclust(mcl.pm, cell_point_size = 0.1)
    
    ## select origin cluster
    pt.pm.mean<- tapply(pt[names(mcl.pm[['clusterid']])], list(mcl.pm[['clusterid']]), mean)
    start.cluster <- names(which.min(pt.pm.mean))
    
    ## construct pseudotime
    ord.pm <- TSCANorder(mcl.pm, startcluster = start.cluster, listbranch = T,orderonly = T)
    # str(ord.pm)
    
    pt.pm <- unlist(sapply(sapply(ord.pm, length), function(i) seq(1, i)))
    names(pt.pm) <- unname(unlist(ord.pm))
    # --- check if these codes are necessary <<<<<<<<<<<<<<<<
    ## plot pseudotime
    
    pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(pt.pm[rownames(pca)]))
    # ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
    #   geom_scattermore() + theme_classic()
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # get candidate branches
    newbranch.pm <- findbranch(mst = mcl.pm$MSTtree, order = ord.pm, origin = start.cluster)
    
    ## compare two MST
    js <- sapply(seq(1, length(newbranch)), function(i){
            id <- which(sapply(paste0(names(ord),','), function(k) grepl(paste0(paste0(newbranch[[i]], collapse = ','),','), k)))[1]
            cells <- ord[[id]]
            b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k) names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
            sapply(seq(1, length(newbranch.pm)), function(j){
              print(j)
              id <- which(sapply(paste0(names(ord.pm),','), function(k) grepl(paste0(paste0(newbranch.pm[[j]], collapse = ','),','), k)))[1]
              cells <- ord.pm[[id]]
              b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k) names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
              js <- length(intersect(b.pm, b.ori))/length(union(b.pm, b.ori))
            })
          })
    oc <- sapply(seq(1, length(newbranch)), function(i){
              id <- which(sapply(paste0(names(ord),','), function(k) grepl(paste0(paste0(newbranch[[i]], collapse = ','),','), k)))[1]
              cells <- ord[[id]]
              b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k) names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
              sapply(seq(1, length(newbranch.pm)), function(j){
                  id <- which(sapply(paste0(names(ord.pm),','), function(k) grepl(paste0(paste0(newbranch.pm[[j]], collapse = ','),','), k)))[1]
                  cells <- ord.pm[[id]]
                  b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k) names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
                  oc <- length(intersect(b.pm, b.ori))/min(length(b.pm), length(b.ori))
             }) 
          })
    corr <- sapply(seq(1, length(newbranch)), function(i){
                id <- which(sapply(paste0(names(ord),','), function(k) grepl(paste0(paste0(newbranch[[i]], collapse = ','),','), k)))[1]
                cells <- ord[[id]]
                b.ori <- intersect(unlist(sapply(newbranch[[i]], function(k) names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
                
                sapply(seq(1, length(newbranch.pm)), function(j){
                    id <- which(sapply(paste0(names(ord.pm),','), function(k) grepl(paste0(paste0(newbranch.pm[[j]], collapse = ','),','), k)))[1]
                    cells <- ord.pm[[id]]
                    b.pm <- intersect(unlist(sapply(newbranch.pm[[j]], function(k) names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
                    ov = intersect(b.ori, b.pm)
                    cor(pt[ov], pt.pm[ov])
                }) 
            })
    corr[is.na(corr)] <- 0
    colnames(corr) <- colnames(oc) <- colnames(js) <- paste0('original', seq(1, length(newbranch)))
    
    ## get js binary to match branches 
    js.binary <- get_binary(js, js.cut)
    corr.score[[pmid]] <- corr * js.binary
    js.melt <- melt(js.binary)
    js.melt <- js.melt[js.melt[,3]!=0,]
    colnames(js.melt) <- c('permutation.branch','original.branch','matched')
    reproduce.js[[pmid]] <- as.character(js.melt[,2])
    
    ## get oc binary to match branches
    oc.binary <- get_binary(oc, oc.cut)
    oc.melt <- melt(oc.binary)
    oc.melt <- oc.melt[oc.melt[,3]!=0,]
    reproduce.oc[[pmid]] <- as.character(oc.melt[,2])
  
    ## samples cell compositions 
    ctcomp <- sapply(js.melt[,2], function(tmp){
      c <- names(clu)[clu %in% newbranch.pm[[tmp]]]
      ctcomp <- rep(0, length(unique(alls)))
      names(ctcomp) <- unique(alls)
      ctcomp[names(table(alls[c]))] <- table(alls[c])
    })
    colnames(ctcomp) <- paste0('origin', js.melt[,2])
    ctcomp <- ctcomp/rowSums(ctcomp)
    
    ctcomp.new <- matrix(0, nrow = length(unique(alls)), ncol = length(newbranch))
    colnames(ctcomp.new) <- paste0('origin', seq(1, length(newbranch)))
    rownames(ctcomp.new) <- unique(alls)
    ctcomp.new[rownames(ctcomp), colnames(ctcomp)] <- ctcomp
    ctcomplist[[pmid]] <- t(ctcomp.new)
  }
  
  reproduce.js <- unlist(reproduce.js)  
  js.perc <- rep(0, length(newbranch))
  js.perc[as.numeric(names(table(reproduce.js)))] <-  table(reproduce.js)/n.permute
  names(js.perc) <- newbranch
  
  reproduce.oc <- unlist(reproduce.oc)  
  oc.perc <- rep(0, length(newbranch))
  oc.perc[as.numeric(names(table(reproduce.oc)))] <-  table(reproduce.oc)/n.permute
  names(oc.perc) <- newbranch
  
  corr.score.m <- do.call(rbind, corr.score)
  corr.score.v <- colSums(corr.score.m)/n.permute
  names(corr.score.v) <- newbranch
  
  sort((js.perc + oc.perc)/2)
  
  detection.rate <- data.frame(detection.rate = (js.perc + oc.perc[names(js.perc)])/2, stringsAsFactors = FALSE)
  sample.cellcomp.mean <- apply(simplify2array(ctcomplist), 1:2, mean)
  sample.cellcomp.sd <- apply(simplify2array(ctcomplist), 1:2, sd)
  rownames(sample.cellcomp.mean) <- newbranch[as.numeric(sub('origin', '', rownames(sample.cellcomp.mean)))]
  rownames(sample.cellcomp.sd) <- newbranch[as.numeric(sub('origin', '', rownames(sample.cellcomp.sd)))]
  
  result <- list(detection.rate = detection.rate, 
                 sample.cellcomp.mean = sample.cellcomp.mean, 
                 sample.cellcomp.sd = sample.cellcomp.sd)
  return(result)
}

result <- evaluate_uncertainty(a, 3)

