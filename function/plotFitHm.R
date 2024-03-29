plotFitHm <- function(testobj, showRowName = FALSE, cellWidthTotal = 250, cellHeightTotal = 400, showCluster = FALSE, colann = NULL, rowann = NULL, annotation_colors = NULL, type = 'TIME', subsampleCell = TRUE, numSubsampleCell=1e3){
  ## cellHeightTotal: when showRowName = TRUE, cellHeightTotal is suggested to be ten times the number of genes (rows).
  ## showCluster: (no implemented yet). if TRUE, "cluster" should be a slot in testobj, and it will be label in the heatmap. If FALSE, no need to pass in "cluster". 
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(ggplot2)
  fit <- testobj$populationFit
  if (subsampleCell){
    if (toupper(type) == 'TIME'){
      id <- round(seq(1, ncol(fit), length.out = numSubsampleCell))
      fit <- fit[, id]
    } else if (toupper(type) == 'VARIABLE'){
      id <- round(seq(1, ncol(fit[[1]]), length.out = numSubsampleCell))
      for (i in 1:length(fit)){
        fit[[i]] <- fit[[i]][, id]
      }
    }
    testobj$pseudotime <- sort(sample(testobj$pseudotime, numSubsampleCell))
    rownames(testobj$cellanno) <- testobj$cellanno[,1]
    testobj$cellanno <- testobj$cellanno[names(testobj$pseudotime), ]
    print('subsample done!')
  } 
    if ('expr.ori' %in% names(testobj)){
      expr <- testobj$expr.ori[, names(testobj$pseudotime)]  
    } else {
      expr <- testobj$expr[, names(testobj$pseudotime)]
    }

  fit.bak = fit
  clu <- testobj$cluster
  if (toupper(type) == 'VARIABLE'){
    if ('DDGType' %in% names(testobj)) DDGType <- testobj$DDGType else DDGType <- getDDGType(testobj) 
    fit.scale <- do.call(cbind, fit)
    fit.scale <- fit.scale[names(testobj$cluster), ]
    fit.scale <- scalematrix(fit.scale)
    colnames(fit.scale) <- paste0(rep(names(fit), each = ncol(fit.scale)/length(fit)), ';cell', seq(1, ncol(fit.scale)))
    res <- data.frame(clu = clu, 
                      cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(1, ncol(fit.scale)/2)], seq(1, ncol(fit.scale)/2))),
                      # changepoint = sapply(names(clu), function(i) which.min(abs(fit.scale[i, seq(1, ncol(fit.scale)/2)]))),
                      changepoint = sapply(names(clu), function(i) {
                        v <- fit.scale[i, seq(1, ncol(fit.scale)/2)]
                        which(v[-length(v)] * v[-1] < 0)[1]
                        }),
                      DDGType = DDGType[names(clu)])
    res <- res[order(res$clu, res$DDGType, res$changepoint, res$cor), ]
  } else {
    fit.scale <- scalematrix(fit)
    dimnames(fit.scale) <- dimnames(fit)
    res <- data.frame(clu = clu, 
                      cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(1, ncol(fit.scale)/2)], seq(1, ncol(fit.scale)/2))),
                      changepoint = sapply(names(clu), function(i) {
                        v <- fit.scale[i, seq(round(ncol(fit.scale)*0.05), round(ncol(fit.scale)*0.95))]
                        which(v[-length(v)] * v[-1] < 0)[1]
                        }))
                      # changepoint = sapply(names(clu), function(i) which.min(abs(fit.scale[i, seq(round(ncol(fit.scale)*0.01), round(ncol(fit.scale)*0.99))]))))
    
    res <- res[order(res$clu, res$changepoint, res$cor), ]
  }
  
  fit.scale <- fit.scale[rownames(res), ]
  # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
  ## ------------------------
  ## plot original expression 
  ## ------------------------
  cellanno <- testobj$cellanno
  
  expr <- expr[, names(testobj$pseudotime)]
  
  if (toupper(type) == 'VARIABLE'){
    tmp <- lapply(names(fit), function(i){
      expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(testobj$design)[testobj$design[, 2] == sub('.*_','', i)], 1]]
    })
    expr.scale <- do.call(cbind, tmp)
  } else if (toupper(type) == 'TIME'){
    expr.scale <- expr
  }
  expr.scale <- scalematrix(expr.scale)
  expr.scale <- expr.scale[rownames(fit.scale), ]  
  # 
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
    if (type == 'variable'){
      colann <- data.frame(
        # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
        pseudotime = testobj$pseudotime[colnames(expr.scale)],
        group = as.character(testobj$design[cellanno[match(colnames(expr.scale), cellanno[,1]),2],2]),
        expression = 'Original',
        stringsAsFactors = F)
      
      col.group = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$group)))
      names(col.group) = unique(colann$group)
    } else if (type == 'time'){
      colann <- data.frame(
        # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
        pseudotime = testobj$pseudotime[colnames(expr.scale)],
        expression = 'Original',
        stringsAsFactors = F)
    }
  }
  rownames(colann) = colnames(expr.scale)
  col.expression = brewer.pal(n = 8, name = "Pastel1")[1:2]
  names(col.expression) = c('Original', 'Model Fitted')
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  
  if (is.null(rowann)){
    if (type == 'variable'){
      rowann = data.frame(
        cluster = as.character(clu),
        DDGType = as.character(DDGType[names(clu)]),
        stringsAsFactors = F)
    } else if (type == 'time'){
      rowann = data.frame(
        cluster = as.character(clu),
        stringsAsFactors = F)
    }
    rownames(rowann) = names(clu)
    rowann <- rowann[rownames(fit.scale), ,drop=F]
    if (length(unique(clu)) < 8){
      col.clu = brewer.pal(8, 'Set1')[1:length(unique(clu))]
    } else {
      col.clu = colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(clu)))
    }
    names(col.clu) = sort(unique(clu))
  }
  
  if (is.null(colann)| is.null(annotation_colors)){
    if (type == 'variable'){
      col.DDGType = brewer.pal(8, 'Dark2')[1:length(unique(res$DDGType))]
      names(col.DDGType) = unique(res$DDGType)
      annotation_colors = list(
        pseudotime = col.pseudotime,
        group = col.group,
        expression = col.expression,
        cluster = col.clu,
        DDGType = col.DDGType)
    } else if (type == 'time'){
      annotation_colors = list(
        pseudotime = col.pseudotime,
        expression = col.expression,
        cluster = col.clu)
    }
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
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal / ncol(expr.scale),
    cellheight = cellHeightTotal / nrow(expr.scale),
    border_color = NA, silent = TRUE)
  plist[[1]] <- p1[[4]]
  
  ## --------------------
  ## plot fitting values
  ## --------------------
  if (type == 'variable'){
    colann.fit <-data.frame(pseudotime = rep(seq(1, ncol(fit.scale)/length(fit)), length(fit)),
                            group = gsub(sub('_.*', '_', names(fit)[1]),'',sub(';.*', '', colnames(fit.scale))), 
                            expression = 'Model Fitted',
                            stringsAsFactors = F)
  } else if (type == 'time'){
    colann.fit <-data.frame(pseudotime = testobj$pseudotime[colnames(fit.scale)],
                            expression = 'Model Fitted',
                            stringsAsFactors = F)
    col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime)))
    names(col.pseudotime) = unique(colann.fit$pseudotime)
    annotation_colors$pseudotime <- col.pseudotime
  }
  
  rownames(colann.fit) = colnames(fit.scale)
  p2 <- pheatmap(
    fit.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann.fit,
    annotation_row = rowann,
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal / ncol(fit.scale),
    cellheight = cellHeightTotal / nrow(fit.scale),
    border_color = NA, silent = TRUE)
  plist[[3]] <- p2[[4]]
  plist[[2]] <- ggplot(data=NULL) + geom_blank() + theme_void()
  # png(paste0('g.png'),width = 4300,height = 3200,res = 300)
  print(grid.arrange(grobs = plist,layout_matrix=matrix(c(1,1,1,1,2,3,3,3,3),nrow=1)))
  # dev.off()
}  





