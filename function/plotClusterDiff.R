plotClusterDiff <- function(testobj, 
                            gene = names(testobj$cluster),
                            cluster = testobj[['cluster']],
                            each = FALSE,
                            sep = '',
                            reverse = F){
  if ('covariateGroupDiff' %in% names(testobj)){
    fit <- testobj$covariateGroupDiff[gene, ,drop=FALSE]
  } else {
    fit <- getCovariateGroupDiff(testobj = testobj, gene = gene, reverse = reverse)
  }
  colnames(fit) <- seq(1, ncol(fit))
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
if (each){
  pd <- melt(fit)
  colnames(pd) <- c('gene', 'pseudotime', 'covariateGroupDiff')
  pd[,1] <- sub(sep, '', pd[,1])
  p <- ggplot(data = pd) + geom_line(aes(x = pseudotime, y = covariateGroupDiff))+
    theme_classic()+
    facet_wrap(~gene)
}  else {
  clu = cluster[gene]
  tmp <- sapply(sort(unique(clu)), function(i){
    m <- colMeans(fit[clu == i, , drop = FALSE])
  })
  colnames(tmp) <- sort(unique(clu))
  pd <- melt(tmp)
  
  colnames(pd) <- c('pseudotime', 'cluster', 'covariateGroupDiff')
  pd$cluster <- factor(pd$cluster)
  p <- ggplot(data = pd) + geom_line(aes(x = pseudotime, y = covariateGroupDiff, color = cluster))+
    theme_classic() + scale_x_continuous(breaks=c(min(pd$pseudotime),max(pd$pseudotime)))
  if (length(unique(pd$cluster)) < 8){
    p <- p + scale_color_brewer(palette = 'Dark2')
  } else {
    p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$cluster))))
  }
}
   print(p)
}


