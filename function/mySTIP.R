mySTIP <- function(fit, gl) {
  ## fit: gene by cell matrix. colnames: pseudotime of the cells. rownames: gene names.
  ## gl: a vector of genes (should be within the rownames) to be highlighted by texts.
  dn <- dimnames(fit)
  fit <- t(apply(fit, 1, scale))
  dimnames(fit) <- dn
  gene <- row.names(fit)
  gl <- intersect(gl,gene)
  zpdirection <- fit[, ncol(fit)/20] < fit[, ncol(fit)/20*19] ## direction (low to high; or high to low)
  zp <- apply(fit, 1, function(sf) {
    ## for GBM myeloid 
    # setdiff(names(which(sapply(1:(length(sf) - 1), function(i)
    #   sf[i] * sf[i + 1] < 0))),colnames(fit)[c(1:(ncol(fit)/20),(ncol(fit)/20*changePointCut):ncol(fit))])
    names(which(sapply(1:(length(sf) - 1), function(i)
      sf[i] * sf[i + 1] < 0)))
  }) ## change point positions
  zp <- sapply(zp, function(i) {
    if (length(i) > 2) {
      i[1:min(2,length(i))]
    } else {
      i
  }})
  zpnum <- sapply(zp, length)
  inczp <- names(which(zpdirection[zpnum == 1]))
  deczp <- names(which(!zpdirection[zpnum == 1]))
  multipoint <- names(zpnum)[zpnum > 1]
  m1 <- names(which(fit[multipoint, 1] > 0)) ## decrease and then increase
  m2 <- names(which(fit[multipoint, 1] < 0)) ## increase and then decrease
  
  geneorder <- NULL
  
  if (length(deczp) > 0) {
    tmp <- unlist(zp[deczp])
    n <- names(tmp)
    tmp <- match(tmp,colnames(fit))
    names(tmp) <- n
    geneorder <-
      c(geneorder, names(sort(tmp, decreasing = F)))
  }
  if (length(m2) > 0) {
    geneorder <-
      c(geneorder, names(sort(sapply(zp[m2], function(i)
        match(i[1],colnames(fit))))))  
  }
  
  if (length(inczp) > 0) {
    tmp <- unlist(zp[inczp])
    n <- names(tmp)
    tmp <- match(tmp,colnames(fit))
    names(tmp) <- n
    geneorder <- c(geneorder, names(sort(tmp)))
  }
  if (length(m1) > 0) {
    geneorder <-
      c(geneorder, names(sort(sapply(zp[m1], function(i)
        match(i[1],colnames(fit))))))
  }
  geneorder <- rev(geneorder)
  plotdata <- fit[geneorder,]
  plotdata <- melt(plotdata)
  colnames(plotdata) <- c("Gene", "Pseudotime", "Expression")
  plotdata$Pseudotime <- match(plotdata$Pseudotime,colnames(fit))
  gl <- gl[order(match(gl, levels(plotdata$Gene)))]
  p1 <-
    ggplot(plotdata, aes(Pseudotime, Gene, fill = Expression)) + 
    geom_tile() + 
    theme_classic() + 
    scale_fill_gradient2(low ="blue", high = "firebrick", midpoint = 0) + 
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),axis.text.y = element_blank(), axis.line.y = element_blank(), plot.margin = margin(2, 0, 2, 0), axis.title.x = element_text(color = 'black')) + scale_x_continuous(expand = c(0.001, 0.001))
  
  yax <- rep("A", length(gene))
  if (length(gl) < length(gene)) {
    yaxglid <-
    round(seq(length(gene)*0.02, length(gene) * 0.98, length.out = length(gl)))
  } else {
    yaxglid <-
    round(seq(1, length(gene), length.out = length(gl)))
  }
  yax[yaxglid] <- gl
  if (length(yaxglid)<length(yax)) yax[setdiff(1:length(yax), yaxglid)] <- setdiff(gene, gl)
  
  p2 <-
    ggplot() + geom_point(data = data.frame(gene = factor(yax,levels=yax), x =1),
                          aes(x = x, y = gene),
                          col = "white") + geom_text(
                            data = data.frame(
                              text = factor(gl,levels=gl),
                              id = factor(gl,levels=gl),
                              x = 1),
                            aes(x = x, y = id, label = text),
                            size = 3
                          ) + geom_segment(
                            data = data.frame(
                              x = 1.5,
                              xend = 2,
                              y = factor(gl,levels=gl),
                              yend = yax[match(gl, geneorder)]
                            ),
                            aes(
                              x = x,
                              y = y,
                              xend = xend,
                              yend = yend
                            ),
                            size = 0.1
                          ) + theme_classic() + coord_cartesian(xlim = c(0.5, 2)) + theme(
                            axis.title = element_text(color = "white"),
                            axis.line = element_line(color = "white"),
                            axis.text = element_text(color = "white"),
                            axis.ticks = element_line(color = "white"),
                            plot.margin = margin(2, 0, 2, 2)
                          ) + scale_x_continuous(expand = c(0, 0))
  gridExtra::grid.arrange(p2, p1, nrow = 1, layout_matrix = cbind(1, 1,1, 2, 2, 2,2,2,2,2))
}



