infer_tree_structure <- function(pca, ct, origin.celltype, number.cluster = NA, plotdir = getwd(), xlab = 'PC1', ylab = 'PC2', max.clunum=50){
  ## ct: dataframe/matrix, first column is cell name, second column is cell type, third column is sample.
  library(igraph)
  alls <- ct[,3]
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
  # clu <- mykmeans(pr, number.cluster = number.cluster, maxclunum = 50, seed = i)$cluster
  clu <- mykmeans(pr, maxclunum = 50, number.cluster = number.cluster)$cluster
  table(clu)
  pd = data.frame(x = pr[,1], y = pr[,2], clu = as.factor(clu[rownames(pr)]))
  mypalette = colorRampPalette(brewer.pal(9,'Set1'))
  pdf(paste0(plotdir, 'cluster.pdf'), width = (0.7*max(clu)), height = (0.5*max(clu)))
  print(ggplot(data = pd, aes(x = x, y = y, color = clu)) +
    geom_scattermore()+
    scale_color_manual(values = mypalette(max(clu)))+
    theme_classic() + xlab(xlab) + ylab(ylab))
  dev.off()
  ## cell type composition in clusters
  pd = cbind(pd, celltype = ct[match(rownames(pd), ct[,1]),2])
  tab <- table(pd[,3:4])
  tab <- tab/rowSums(tab)
  pd <- melt(tab)
  pd$clu <- factor(as.character(pd$clu), levels = seq(1,max(pd$clu)))
  pdf(paste0(plotdir, 'celltype_composition_for_cluster.pdf'), width = 9, height = 5)
  print(ggplot(data = pd) +
    geom_bar(aes(x = clu, y = value, fill = celltype), stat = 'identity', position = 'dodge') +
    theme_classic() +
    ylab('Celltype Proportion') +
    scale_fill_manual(values = mypalette(length(unique(pd$celltype)))))
  dev.off()
  ### mclust
  mcl <- exprmclust(t(pr),cluster=clu,reduce=F)
  # mcl <- exprmclust(t(pr), reduce = F)
  # pdf(paste0(plotdir, 'mcl.pdf'), width=(0.8*max(clu)),height=(0.6*max(clu)))
  # print(plotmclust(mcl, cell_point_size = 0.1, x.lab = xlab, y.lab = ylab))
  # dev.off()

  # str(mcl)
  # 
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
  pd = data.frame(pc1 = pca[,1], pc2 = pca[,2], time = as.numeric(pt[rownames(pca)]))
  library(scattermore)
  library(RColorBrewer)
  pdf(paste0(plotdir, 'pseudotime.pdf'), width = 5, height = 4.2)  
  print(ggplot(data = pd, aes(x = pc1, y = pc2, color = time)) +
    geom_scattermore() +
    scale_color_gradient(low = 'yellow', high = 'blue')+
    xlab(xlab) + ylab(ylab) +
    theme_classic())
  dev.off()
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



