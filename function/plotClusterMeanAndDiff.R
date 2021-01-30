plotClusterMeanAndDiff <- function(testobj, 
                            cluster = testobj[['cluster']],
                            free.scale = TRUE){
  ## only works for Covariate Test. 
  if ('populationFit' %in% names(testobj)){
    fit <- testobj$populationFit
  } else {
    print("The object testobj should contain populationFit!")
  }
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  library(gridExtra)
  a <- ifelse(free.scale, 'free', 'fixed') 
  int <- intersect(rownames(fit[[1]]), names(cluster))
  clu <- cluster[int]
  pd <- lapply(1:length(fit), function(i){
    mat <- fit[[i]][int, ]
    tmp <- sapply(sort(unique(clu)), function(i){
      colMeans(mat[clu == i, , drop = FALSE])
    })
    tmp <- melt(tmp)
    colnames(tmp) <- c('pseudotime', 'cluster', 'populationFitClusterMean')
    tmp <- data.frame(tmp, type = names(fit[i]))
  })
  pd <- do.call(rbind, pd)
  pd$cluster <- factor(pd$cluster)
  p1 <- ggplot(data = pd) + geom_line(aes(x = pseudotime, y = populationFitClusterMean, color = type), size = 1)+
    theme_classic() + 
    facet_wrap(~cluster, nrow = length(unique(pd$cluster)), scales = a)+
    theme(legend.position = 'none') 
  if (length(unique(pd$type)) < 8){
    p1 <- p1 + scale_color_brewer(palette = 'Dark2')
  } else {
    p1 <- p1 + scale_color_manual(values = colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$type))))
  }
   
   
  if ('covariateGroupDiff' %in% names(testobj)){
    fit <- testobj$covariateGroupDiff
  } else {
    fit <- getCovariateGroupDiff(testobj = testobj, gene = int)
  }
  colnames(fit) <- seq(1, ncol(fit))
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)

  tmp <- sapply(sort(unique(clu)), function(i){
    m <- colMeans(fit[clu == i, , drop = FALSE])
  })
  colnames(tmp) <- sort(unique(clu))
  pd2 <- melt(tmp)
  
  colnames(pd2) <- c('pseudotime', 'cluster', 'covariateGroupDiff')
  pd2$cluster <- factor(pd2$cluster)
  p2<- ggplot(data = pd2) + geom_line(aes(x = pseudotime, y = covariateGroupDiff), size = 1)+
    theme_classic() +
    facet_wrap(~cluster, nrow = length(unique(pd2$cluster)), scales = a)
  if (length(unique(pd$cluster)) < 8){
    p2 <- p2 + scale_color_brewer(palette = 'Set1')
  } else {
    p2 <- p2 + scale_color_manual(values = colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$cluster))))
  }
  grid.arrange(p1,p2,ncol=2)
}


