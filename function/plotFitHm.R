plotFitHm <- function(testobj, showRowName = FALSE, cellWidthTotal = 250, cellHeightTotal = 450, showCluster = FALSE, colann = NULL, rowann = NULL){
  ## cellHeightTotal: when showRowName = TRUE, cellHeightTotal is suggested to be ten times the number of genes (rows).
  ## showCluster: (no implemented yet). if TRUE, "cluster" should be a slot in testobj, and it will be label in the heatmap. If FALSE, no need to pass in "cluster". 
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(ggplot2)
  fit <- testobj$populationFit
  fit.bak = fit
  clu <- testobj$cluster
  if (DEGType %in% names(testobj)) 
    DEGType <- testobj$DEGType
  else 
    DEGType <- getDEGType(testobj) 
  
  fit.scale <- do.call(cbind, fit)
  fit.scale <- fit.scale[names(testobj$cluster), ]
  fit.scale <- scalematrix(fit.scale)
  str(fit.scale)
  colnames(fit.scale) <- seq(1, ncol(fit.scale))
  
  res <- data.frame(clu = clu, 
                    cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(1, ncol(fit.scale)/2)], seq(1, ncol(fit.scale)/2))))
  fit.scale <- fit.scale[rownames(res)[order(res$clu, res$cor)], ]
  # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
  ## ------------------------
  ## plot original expression 
  ## ------------------------
  cellanno <- testobj$cellanno
  expr = testobj$expr.ori
  expr <- expr[, names(testobj$pseudotime)]
  expr.scale <-
    cbind(expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(testobj$design[testobj$design[, 2] == sub('.*_','',names(fit)[1]), ]), 1]],
          expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(testobj$design[testobj$design[, 2] == sub('.*_','',names(fit)[2]), ]), 1]])
  
  ## plot ------------------------
  expr.scale[expr.scale > quantile(as.vector(expr.scale), 0.98)] <-
    quantile(as.vector(expr.scale), 0.98)
  expr.scale[expr.scale < quantile(as.vector(expr.scale), 0.08)] <-
    quantile(as.vector(expr.scale), 0.02)
  fit.scale[fit.scale > quantile(as.vector(fit.scale), 0.98)] <-
    quantile(as.vector(fit.scale), 0.98)
  fit.scale[fit.scale < quantile(as.vector(fit.scale), 0.02)] <-
    quantile(as.vector(fit.scale), 0.02)
  ### annotate rows and columns
  
  
  if (is.null(colann)){
    colann <- data.frame(
        # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
        pseudotime = testobj$pseudotime[colnames(expr.scale)],
        expression = 'Original',
        stringsAsFactors = F)
    rownames(colann) = colnames(expr.scale)
    col.expression = brewer.pal(n = 8, name = "Pastel1")[1:2]
    names(col.expression) = c('Original', 'Model Fitted')
    col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
    names(col.pseudotime) = unique(colann$pseudotime)
  }
    
  if (is.null(rowann)){
    rowann = data.frame(
      cluster = as.character(clu),
      DEGType = as.character(DEGType[names(clu)]),
      stringsAsFactors = F)
    rownames(rowann) = names(clu)
    rowann <- rowann[rownames(fit.scale), ,drop=F]
    if (length(unique(clu)) < 8){
      col.clu = brewer.pal(8, 'Set1')[1:length(unique(clu))]
    } else {
      col.clu = colorRampPalette(brewer.pal(8, 'Set1'))[1:length(unique(clu))]
    }
    names(col.clu) = unique(clu)
    col.DEGType = brewer.pal(8, 'Dark2')[1:length(unique(DEGType))]
    names(col.DEGType) = unique(DEGType)
    
  }
    
  #### save png
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  plist <- list()
  p1 <- pheatmap(
    expr.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann,
    annotation_row = rowann,
    annotation_colors = list(
      pseudotime = col.pseudotime,
      expression = col.expression,
      cluster = col.clu,
      DEGType = col.DEGType),
    cellwidth = cellWidthTotal / ncol(expr.scale),
    cellheight = cellHeightTotal / nrow(expr.scale),
    border_color = NA, silent = TRUE)
  plist[[1]] <- p1[[4]]
  
  ## --------------------
  ## plot fitting values
  ## --------------------
  colann2 <-data.frame(pseudotime = rep(seq(1, ncol(fit.scale)/length(fit)), length(fit)),
               expression = 'Model Fitted',
               stringsAsFactors = F)
  rownames(colann2) = colnames(fit.scale)
  p2 <- pheatmap(
    fit.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann2,
    annotation_row = rowann,
    annotation_colors = list(
      expression = col.expression,
      pseudotime = col.pseudotime,
      cluster = col.clu,
      DEGType = col.DEGType),
      cellwidth = cellWidthTotal / ncol(fit.scale),
      cellheight = cellHeightTotal / nrow(fit.scale),
      border_color = NA, silent = TRUE)
  plist[[3]] <- p2[[4]]
  plist[[2]] <- ggplot(data=NULL) + geom_blank() + theme_void()
  
  # png(paste0('g.png'),width = 4300,height = 3200,res = 300)
  g <- grid.arrange(grobs = plist,layout_matrix=matrix(c(1,1,1,1,2,3,3,3,3),nrow=1))
  # dev.off()
  return(g)
}  
  




