plotClusterMean <- function(testobj, 
                            cluster,
                            type = 'time',
                            facet = FALSE, 
                            facet_scales = 'free',
                            facet_nrow = 3){
  ## type: "time" (default) or "variable.
  ## facet: plot each cluster indiviidually. Lnly works when type == 'time'. When type == 'variable', the plot will be facet anyway for each cluster.
  if ('populationFit' %in% names(testobj)){
    fit <- testobj$populationFit
  } else {
    print("The object testobj should contain populationFit!")
  }
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  if (toupper(type) == 'VARIABLE'){
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
      theme(axis.text.x = element_blank()) +
      facet_wrap(~cluster, scales = facet_scales, nrow = facet_nrow) +
      xlab('Pseudotime') +
      ylab('Population fitting cluster mean')
    if (length(unique(pd$type)) < 8){
      p <- p + scale_color_brewer(palette = 'Dark2')
    } else {
      p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$type))))
    }
  } else if (toupper(type) == 'TIME'){
    int <- intersect(rownames(fit), names(cluster))
    clu <- cluster[int]
    mat <- fit[int, ]
    tmp <- sapply(sort(unique(clu)), function(i){
      colMeans(mat[clu == i, , drop = FALSE])
    })
    pd <- melt(tmp)
    colnames(pd) <- c('pseudotime', 'cluster', 'populationFitClusterMean')
    pd$pseudotime <- as.numeric(pd$pseudotime)
    pd$cluster <- factor(pd$cluster, levels = sort(unique(pd$cluster)))
    p <- ggplot(data = pd, aes(x = pseudotime, y = populationFitClusterMean, group = cluster, color = cluster))+
      geom_smooth(size = 1) +
      theme_classic() +
      theme(axis.text.x = element_blank()) +
      xlab('Pseudotime') +
      ylab('Population fitting cluster mean')
    if (facet){
      p <- p + facet_wrap(~cluster, scales = facet_scales, nrow = facet_nrow)
    }
    if (length(unique(pd$cluster)) < 8){
      p <- p + scale_color_brewer(palette = 'Set1')
    } else {
      p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(9,'Set1'))(length(unique(pd$cluster))))
    }
  }
    
   print(p)
}




