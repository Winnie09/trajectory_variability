plotDiffFitHm <- function(testobj, showRowName = FALSE, cellWidthTotal = 250, cellHeightTotal = 400, showCluster = FALSE, colann = NULL, rowann = NULL, annotation_colors = NULL, type = 'time', subsampleCell = TRUE, numSubsampleCell=1e3){
  ## cellHeightTotal: when showRowName = TRUE, cellHeightTotal is suggested to be ten times the number of genes (rows).
  ## showCluster: (no implemented yet). if TRUE, "cluster" should be a slot in testobj, and it will be label in the heatmap. If FALSE, no need to pass in "cluster". 
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(ggplot2)
  fit <- testobj$populationFit
  
  if (subsampleCell){
    if (type == 'time'){
      id <- round(seq(1, ncol(fit), length.out = numSubsampleCell))
      fit <- fit[, id]
    } else if (type == 'variable'){
      id <- round(seq(1, ncol(fit[[1]]), length.out = numSubsampleCell))
      for (i in 1:length(fit)){
        fit[[i]] <- fit[[i]][, id]
      }
      FitDiff.scale <- scalematrix(testobj$covariateGroupDiff[,id]) ## add FitDiff.scale
      colnames(FitDiff.scale) <- paste0('FitDiff:cell', seq(1, ncol(FitDiff.scale)))
      testobj$pseudotime <- sort(sample(testobj$pseudotime, numSubsampleCell))
      rownames(testobj$cellanno) <- testobj$cellanno[,1]
      testobj$cellanno <- testobj$cellanno[names(testobj$pseudotime), ]
      testobj$expr.ori <- testobj$expr.ori[, names(testobj$pseudotime)]
      
    }
  }
  print('subsample done!')
  fit.bak = fit
  clu <- testobj$cluster
  if (type == 'variable'){
    if ('DEGType' %in% names(testobj)) 
      DEGType <- testobj$DEGType
    else 
      DEGType <- getDEGType(testobj) 
    fit.scale <- do.call(cbind, fit)
    fit.scale <- fit.scale[names(testobj$cluster), ]
    fit.scale <- scalematrix(fit.scale)
    colnames(fit.scale) <- paste0(rep(names(fit), each = ncol(fit.scale)/length(fit)), ';cell', seq(1, ncol(fit.scale)))
    
  } else {
    fit.scale <- scalematrix(fit)
    dimnames(fit.scale) <- dimnames(fit)
    
  }
  # res <- data.frame(clu = clu, 
  #                   cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(1, ncol(fit.scale)/2)], seq(1, ncol(fit.scale)/2))),
  #                   changepoint = sapply(names(clu), function(i) which.min(abs(fit.scale[i, seq(1, ncol(fit.scale)/2)]))),
  #                   DEGType = DEGType[names(clu)])
  res <- data.frame(clu = clu, 
                    cor = sapply(names(clu), function(i) cor(FitDiff.scale[i, seq(1, ncol(FitDiff.scale))], seq(1, ncol(FitDiff.scale)))),
                    changepoint = sapply(names(clu), function(i) which.min(abs(FitDiff.scale[i, seq(1, ncol(FitDiff.scale))]))),
                    DEGType = DEGType[names(clu)])
  
  res <- res[order(res$clu, res$DEGType, res$changepoint, res$cor), ]
  fit.scale <- fit.scale[rownames(res), ]
  FitDiff.scale <- FitDiff.scale[rownames(res), ]
  # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
  ## ------------------------
  ## plot original expression 
  ## ------------------------
  cellanno <- testobj$cellanno
  expr = testobj$expr.ori
  expr <- expr[, names(testobj$pseudotime)]
  
  if (type == 'variable'){
    tmp <- lapply(names(fit), function(i){
      expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(testobj$design)[testobj$design[, 2] == sub('.*_','', i)], 1]]
    })
    expr.scale <- do.call(cbind, tmp)
  } else if (type == 'time'){
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
  FitDiff.scale[FitDiff.scale > quantile(as.vector(FitDiff.scale), 0.98)] <-
    quantile(as.vector(FitDiff.scale), 0.98)
  FitDiff.scale[FitDiff.scale < quantile(as.vector(FitDiff.scale), 0.08)] <-
    quantile(as.vector(FitDiff.scale), 0.02)
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
      
      col.group = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$group))+1)
      names(col.group) = c('NA', unique(colann$group))
    } else if (type == 'time'){
      colann <- data.frame(
        # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
        pseudotime = testobj$pseudotime[colnames(expr.scale)],
        expression = 'Original',
        stringsAsFactors = F)
    }
  }
  rownames(colann) = colnames(expr.scale)
  col.expression = brewer.pal(n = 8, name = "Pastel2")[1:3]
  names(col.expression) = c('Original', 'ModelFitted', 'ModeledGroupDiff')
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  
  if (is.null(rowann)){
    if (type == 'variable'){
      rowann = data.frame(
        cluster = as.character(clu),
        DEGType = as.character(DEGType[names(clu)]),
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
      col.clu = colorRampPalette(brewer.pal(8, 'Set1'))[1:length(unique(clu))]
    }
    names(col.clu) = unique(clu)
  }
  
  if (is.null(colann)| is.null(annotation_colors)){
    if (type == 'variable'){
      col.DEGType = brewer.pal(8, 'Set3')[1:length(unique(res$DEGType))]
      names(col.DEGType) = unique(res$DEGType)
      annotation_colors = list(
        pseudotime = col.pseudotime,
        group = col.group,
        expression = col.expression,
        cluster = col.clu,
        DEGType = col.DEGType)
    } else if (type == 'time'){
      annotation_colors = list(
        pseudotime = col.pseudotime,
        expression = col.expression,
        cluster = col.clu,
        gs = col.gs,
        limmaPb = col.limmaPb)
    }
  }
  
  
  col.gs <- c('pink', 'skyblue')
  names(col.gs) <- c('No', 'Yes')
  col.limmaPb <- c('pink', 'skyblue')
  names(col.limmaPb) <- c('Diff', 'nonDiff')
  annotation_colors[['gs']] <- col.gs
  annotation_colors[['limmaPb']] <- col.limmaPb
  
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
    border_color = NA) ## , silent = TRUE
  plist[[1]] <- p1[[4]] 
  
  ## --------------------
  ## plot fitting values
  ## --------------------
  if (type == 'variable'){
    colann.fit1 <-data.frame(pseudotime = rep(testobj$pseudotime, length(fit)),
                             group = gsub(sub('_.*', '_', names(fit)[1]),'',sub(';.*', '', colnames(fit.scale))), 
                             expression = 'ModelFitted',
                             stringsAsFactors = F)
    colann.fit2 <-data.frame(pseudotime = testobj$pseudotime,
                             group = 'NA', 
                             expression = 'ModeledGroupDiff',
                             stringsAsFactors = F)
    colann.fit <- rbind(colann.fit1, colann.fit2)
    
  } else if (type == 'time'){
    colann.fit <-data.frame(pseudotime = testobj$pseudotime[colnames(fit.scale)],
                            expression = 'ModelFitted',
                            stringsAsFactors = F)
  }
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime)))
  names(col.pseudotime) = unique(colann.fit$pseudotime)
  annotation_colors$pseudotime <- col.pseudotime
  # col.group <- colorRampPalette(brewer.pal(n=8,'Accent'))(length(unique(colann.fit$group)))
  # names(col.group) = unique(colann.fit$group)
  # annotation_colors$group <-   col.group
  
  
  fit.scale <- cbind(fit.scale, FitDiff.scale)  ## cbind FitDiff !!!
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
    cellwidth = cellWidthTotal*1.23 / ncol(fit.scale),
    cellheight = cellHeightTotal / nrow(fit.scale),
    border_color = NA) ## , silent = TRUE
  plist[[3]] <- p2[[4]] 
  plist[[2]] <- ggplot(data=NULL) + geom_blank() + theme_void()
  
  # png(paste0('g.png'),width = 4300,height = 3200,res = 300)
  print(grid.arrange(grobs = plist,layout_matrix=matrix(c(1,1,1,1,2,3,3,3,3),nrow=1)))
  # dev.off()
  
}  





