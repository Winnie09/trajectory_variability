mySTIP2 <- function(fit, gl, plot = TRUE, ReturnGeneOrder = FALSE) {
  ## this function will show all gene lables of gl. Heatmap is similar to STIP, but only contains gl genes.
  ## fit: gene by cell matrix. colnames: pseudotime of the cells. rownames: gene names.
  dn <- dimnames(fit)
  fit <- t(apply(fit, 1, scale))
  dimnames(fit) <- dn
  gene <- row.names(fit)
  gl <- intersect(gl,gene)
  
  zpdirection <- fit[, 1] < fit[, ncol(fit)]
  
  zp <- apply(fit, 1, function(sf) {
    names(which(sapply(1:(length(sf) - 1), function(i)
      sf[i] * sf[i + 1] < 0)))
  })
  zpnum <- sapply(zp, length)
  inczp <- names(which(zpdirection[zpnum == 1]))
  deczp <- names(which(!zpdirection[zpnum == 1]))
  multipoint <- names(zpnum)[zpnum > 1]
  m1 <- names(which(fit[multipoint, 1] > 0))
  m2 <- names(which(fit[multipoint, 1] < 0))
  
  geneorder <- NULL
  
  if (length(deczp) > 0) {
    tmp <- unlist(zp[deczp])
    n <- names(tmp)
    tmp <- match(tmp,colnames(fit))
    names(tmp) <- n
    geneorder <-
      c(geneorder, names(sort(tmp, decreasing = F)))
  }
  if (length(m2)>0)
  geneorder <-
    c(geneorder, names(sort(sapply(zp[m2], function(i)
      match(i[1],colnames(fit))))))
  if (length(inczp) > 0) {
    tmp <- unlist(zp[inczp])
    n <- names(tmp)
    tmp <- match(tmp,colnames(fit))
    names(tmp) <- n
    geneorder <- c(geneorder, names(sort(tmp)))
  }
  if (length(m1)>0) geneorder <-
    c(geneorder, names(sort(sapply(zp[m1], function(i)
      match(i[1],colnames(fit))))))
  geneorder <- rev(geneorder)
  if (ReturnGeneOrder) {
    return(geneorder)
  } else {
    plotdata <- fit[geneorder,]
    plotdata <- melt(plotdata)
    colnames(plotdata) <- c("Gene", "Pseudotime", "Expression")
    plotdata$Pseudotime <- match(plotdata$Pseudotime,colnames(fit))
    gl <- gl[order(match(gl, levels(plotdata$Gene)))]
    p <-
      ggplot(plotdata, aes(Pseudotime, Gene, fill = Expression)) + geom_tile() + theme_classic() + scale_fill_gradient2(low ="blue",high = "red",midpoint = 0) + theme(
                                                                                                                          axis.title.y = element_blank(),
                                                                                                                          axis.ticks.y = element_blank(),
                                                                                                                          axis.line.y = element_blank(),
                                                                                                                          plot.margin = margin(2, 0, 2, 0)
                                                                                                                        ) + scale_x_continuous(expand = c(0.001, 0.001))
    return(p)
  }
}



