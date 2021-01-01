plotClusterMean <- function(testobj, 
                            cluster){
  ## only works for Covariate Test. 
  if ('populationFit' %in% names(testobj)){
    fit <- testobj$populationFit
  } else {
    print("The object testobj should contain populationFit!")
  }
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
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
  p <- ggplot(data = pd) + geom_line(aes(x = pseudotime, y = populationFitClusterMean, color = type), size = 1)+
    theme_classic() + 
    facet_wrap(~cluster)
  if (length(unique(pd$type)) < 8){
    p <- p + scale_color_brewer(palette = 'Dark2')
  } else {
    p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$type))))
  }
   print(p)
}



