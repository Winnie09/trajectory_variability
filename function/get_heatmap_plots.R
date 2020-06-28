get_heatmap_plots <- function(testptObj, Mat, Pseudotime, Cellanno, Design, Max.gene = NULL){
  library(ggplot2)
  library(ggdendro)
  library(reshape2)
  if (ncol(Pseudotime) > 1) {
    Pseudotime <- data.frame(Cell = Pseudotime[,1], Pseudotime = as.numeric(Pseudotime[,2]), stringsAsFactors = FALSE)
  Pseudotime <- Pseudotime[order(Pseudotime[,2]), ]
  psn <- Pseudotime[,2]
  names(psn) <- Pseudotime[,1]
  Pseudotime <- psn
} else {
  print('Note: pseudotime should be a dataframe containing 1st column cell and 2nd column pseudotime.')
} 
  Order <- data.frame(Cell = names(Pseudotime), Pseudotime = Pseudotime, stringsAsFactors = FALSE)
  colnames(Cellanno) <- c('Cell', 'Sample')
  Mat <- Mat[, Cellanno[,1], drop=F]
  knotnum <- testptObj$knotnum
  knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
  
  # Find significant genes
  id = intersect(rownames(Mat), names(res$fdr))
  stat <- data.frame(fdr = res$fdr[id], foldchange = res$foldchange[id], stringsAsFactors = FALSE)
  stat <- stat[order(stat$fdr, -(stat$foldchange)), ]
  stat <- stat[stat$fdr <0.05, , drop = FALSE]
  Gene <- rownames(stat)
  if (!is.null(Max.gene)){
    Gene <- Gene[1:Max.gene]
  }
  # prepare the matrix for heatmap, 1000 * num.gene, values are association levels
  expr.time.cor <- sapply(Gene, function(g){
      pd <- data.frame(Expr = Mat[g, ], Cell = Cellanno[,1], Sample = Cellanno[,2], Variable = Design[match(Cellanno[,2], rownames(Design)), 1])
      pd = cbind(pd, Pseudotime = Order[match(pd$Cell, Order$Cell),'Pseudotime'], g = g)
      linedlist <- lapply(unique(pd$Sample), function(p){
        tMat = Mat[g, which(Cellanno[,2]==p),drop=F]
        trainX = Order[match(colnames(tMat), Order$Cell),'Pseudotime']      ### use time 
        pred <- get_spline_fit(tMat, trainX=trainX, fit.min=min(Order$Pseudotime), fit.max=max(Order$Pseudotime), num.base = knotnum[g])
        tmpdf <- data.frame(Expr=pred[1,], Pseudotime=seq(min(Order$Pseudotime),max(Order$Pseudotime),length.out = 1000), Sample=p, Variable = Design[rownames(Design) == p, 1], g = g) ### use here to calculate cor !!!!!!!
      })
      ld = do.call(rbind, linedlist)
      c <- sapply(unique(ld$Pseudotime), function(t){
        tmp = ld[ld$Pseudotime == t, ]
        cor(tmp$Expr, tmp$Variable)
      })
  })  

  # Scale each measurement (independently) to have a mean of 0 and variance of 1
  otter.scaled <- t(scale(expr.time.cor))
  
  # Run clustering
  otter.matrix <- as.matrix(otter.scaled)
  
  otter.dendro <- as.dendrogram(hclust(d = dist(x = otter.matrix)))
  
  # Create dendrogram plot
  dendro.plot <- ggdendrogram(data = otter.dendro, rotate = TRUE) + 
    theme(axis.text.y = element_text(size = 6))
  
  # Heatmap
  
  # Data wrangling
  otter.long <- melt(otter.scaled)
  # Extract the order of the tips in the dendrogram
  otter.order <- order.dendrogram(otter.dendro)
  # Order the levels according to their position in the cluster
  otter.long$Var1 <- factor(x = otter.long$Var1,
                                 levels = rownames(otter.scaled)[otter.order], 
                                 ordered = TRUE)
  # Create heatmap plot
  heatmap.plot <- ggplot(data = otter.long, aes(x = Var2, y = Var1)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "top") + 
    xlab('Pseudotime') +
    ylab('Gene') +
    labs(fill = 'Association Level')
  return(list(heatmap.plot = heatmap.plot, dendro.plot = dendro.plot))
}

