# testobj <- Res
# showRowName = F
# cellWidthTotal = 200
# cellHeightTotal = 300
# # cellWidthTotal = 300
# # cellHeightTotal = length(Res$cluster) * 10
# sep = ':.*'
# subsampleCell = T
# break.0 = FALSE
# showCluster = FALSE
# 
# colann = NULL
# rowann = NULL
# annotation_colors = NULL
# subsampleCell = TRUE
# numSubsampleCell=1e3
# sep = ':.*'
# break.0 = TRUE

plotDiffFitHm4 <- function(testobj, showRowName = FALSE, cellWidthTotal = 250, cellHeightTotal = 400, showCluster = FALSE, colann = NULL, rowann = NULL, annotation_colors = NULL, subsampleCell = TRUE, numSubsampleCell=1e3, sep = NA, break.0 = TRUE){
  ## cellHeightTotal: when showRowName = TRUE, cellHeightTotal is suggested to be ten times the number of genes (rows).
  ## showCluster: (no implemented yet). if TRUE, "cluster" should be a slot in testobj, and it will be label in the heatmap. If FALSE, no need to pass in "cluster". 
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(ggplot2)
  fit <- testobj$populationFit
  if ('DDGType' %in% names(Res)) {
    DDGType <- Res$DDGType
  } else {
    DDGType <- getDDGType(Res)
  }
  DDGType <- DDGType[rownames(testobj$covariateGroupDiff)]
  
  if (subsampleCell){
    id <- round(seq(1, ncol(fit[[1]]), length.out = numSubsampleCell))
    for (i in 1:length(fit)){
      fit[[i]] <- fit[[i]][, id]
    }
    if (sum(DDGType == 'meanSig',na.rm=T) > 0){
      meanid <- which(DDGType == 'meanSig')
      ## <<<<<<<<<<<<<< scale group difference by absmax
      # FitDiff.scale1 <- scalematrix(testobj$covariateGroupDiff[-meanid,id,drop=F]) ## add FitDiff.scale
      max <- apply(abs(testobj$covariateGroupDiff[-meanid,id,drop=F]), 1, max)
      FitDiff.scale1 <- testobj$covariateGroupDiff[-meanid,id,drop=F]/max
      ## >>>>>>>>>>>>>>>
      FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid,id,drop=F])
      FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff),, drop=F]
    } else {
      ## <<<<<<<<<<<<<< scale group difference by absmax
      # FitDiff.scale <- scalematrix(testobj$covariateGroupDiff[,id,drop=F]) ## add FitDiff.scale
      max <- apply(abs(testobj$covariateGroupDiff[,id,drop=F]), 1, max)
      FitDiff.scale <- testobj$covariateGroupDiff[,id,drop=F]/max
      ## >>>>>>>>>>>>>>
    }
    FitDiff.sd <- scalematrix(testobj$covariateGroupDiff[,id,drop=F]) ## add FitDiff.scale  
    colnames(FitDiff.scale) <- paste0('FitDiff:cell', seq(1, ncol(FitDiff.scale)))
    testobj$pseudotime <- sort(sample(testobj$pseudotime, numSubsampleCell))
    print('subsample done!')
  } else {
    if (sum(DDGType == 'meanSig') > 0){
      meanid <- which(DDGType == 'meanSig')
      ## <<<<<<<<<<<<<< scale group difference by absmax
      # FitDiff.scale1 <- scalematrix(testobj$covariateGroupDiff[-meanid,,drop=F]) ## add FitDiff.scale
      max <- apply(abs(testobj$covariateGroupDiff[-meanid,,drop=F]), 1, max)
      FitDiff.scale1 <- testobj$covariateGroupDiff[-meanid,,drop=F]/max
      ## >>>>>>>>>>>>>>
      FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid,,drop=F])
      FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff),, drop=F]
    } else {
      ## <<<<<<<<<<<<<< scale group difference by absmax
      # FitDiff.scale <- scalematrix(testobj$covariateGroupDiff)
      max <- apply(abs(testobj$covariateGroupDiff), 1, max)
      FitDiff.scale <- testobj$covariateGroupDiff/max
      ## >>>>>>>>>>>>>>
      
    }
    colnames(FitDiff.scale) <- paste0('FitDiff:cell', seq(1, ncol(FitDiff.scale)))
    FitDiff.sd <- scalematrix(testobj$covariateGroupDiff) ## add FitDiff.scale  
  }
  oridata <- testobj$covariateGroupDiff
  
  fit.bak = fit
  clu <- testobj$cluster
  
  rownames(testobj$cellanno) <- testobj$cellanno[,1]
  testobj$cellanno <- testobj$cellanno[names(testobj$pseudotime), ]
  if ('expr.ori' %in% names(testobj)){
    testobj$expr <- testobj$expr.ori[, names(testobj$pseudotime)]
  } else {
    testobj$expr <- testobj$expr[, names(testobj$pseudotime)]
  }
  
  if ('DDGType' %in% names(testobj)) {
    DDGType <- testobj$DDGType
  } else {
    DDGType <- getDDGType(testobj) 
  }
  
  fit.scale <- do.call(cbind, fit)
  fit.scale <- fit.scale[names(testobj$cluster), ]
  fit.scale <- scalematrix(fit.scale)
  colnames(fit.scale) <- paste0(rep(names(fit), each = ncol(fit.scale)/length(fit)), ';cell', seq(1, ncol(fit.scale)))
  
  res <- data.frame(clu = clu, 
                    cor = sapply(names(clu), function(i) cor(FitDiff.scale[i, seq(1, ncol(FitDiff.scale))], seq(1, ncol(FitDiff.scale)))),
                    # changepoint = sapply(names(clu), function(i) which.min(abs(FitDiff.sd[i, seq(1, ncol(FitDiff.scale))]))),
                    changepoint = sapply(names(clu), function(i){which(FitDiff.sd[i, -ncol(FitDiff.sd)] * FitDiff.sd[i, -1] < 0)[1]}),
                    DDGType = DDGType[names(clu)])
  
  #res <- res[order(res$clu, res$DDGType, res$changepoint, res$cor), ]
  res1 <- res[res$DDGType!='meanSig',]
  res2 <- res[res$DDGType=='meanSig',]
  
  o1 <- rownames(res1)[order(match(res1$DDGType,c('trendSig','bothSig','other')), res1$cor > 0,res1$changepoint)]
  o2 <- rownames(res2)[order(res2$clu)]
  res <- res[c(o1,o2), ]
  fit.scale <- fit.scale[rownames(res), ]
  FitDiff.scale <- FitDiff.scale[rownames(res), ]
  # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
  ## ------------------------
  ## plot original expression 
  ## ------------------------
  cellanno <- testobj$cellanno
  expr = testobj$expr
  expr <- expr[, names(testobj$pseudotime)]
  
  
  tmp <- lapply(names(fit), function(i){
    expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(testobj$design)[testobj$design[, 2] == sub('.*_','', i)], 1]]
  })
  expr.scale <- do.call(cbind, tmp)
  expr.scale <- scalematrix(expr.scale)
  expr.scale <- expr.scale[rownames(fit.scale), ]  
  
  ## plot data 
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
    colann <- data.frame(
      # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
      pseudotime = testobj$pseudotime[colnames(expr.scale)],
      group = as.character(testobj$design[cellanno[match(colnames(expr.scale), cellanno[,1]),2],2]),
      expression = 'Original',
      stringsAsFactors = F)
    
    col.group = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$group))+1)
    names(col.group) = c('NA', unique(colann$group))
  }
  rownames(colann) = colnames(expr.scale)
  col.expression = brewer.pal(n = 8, name = "Pastel2")[1:3]
  names(col.expression) = c('Original', 'ModelFitted', 'ModeledGroupDiff')
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  
  if (is.null(rowann)){
    rowann = data.frame(
      cluster = as.character(clu),
      DDGType = as.character(DDGType[names(clu)]),
      stringsAsFactors = F)
    rownames(rowann) = names(clu)
  }
  rowann <- rowann[rownames(fit.scale), ,drop=F]
  rowann[,'DDGType'] <- factor(as.character(rowann[,'DDGType']), levels = c('trendSig','meanSig','bothSig','nonDDG','other'))
  
  if (length(unique(clu)) < 8){
    col.clu = brewer.pal(8, 'Set1')[1:length(unique(clu))]
  } else {
    col.clu = colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(clu)))
    
  }
  names(col.clu) = unique(clu)
  
  if (is.null(colann)| is.null(annotation_colors)){
    
    col.meanDiff = c('blue','red')
    names(col.meanDiff) <- c('Positive','Negative')
    col.DDGType = brewer.pal(8, 'Set3')[1:5]
    names(col.DDGType) = c('trendSig','meanSig','bothSig','nonDDG','other')
    annotation_colors = list(
      pseudotime = col.pseudotime,
      group = col.group,
      expression = col.expression,
      cluster = col.clu,
      DDGType = col.DDGType,
      meanDiff = col.meanDiff)
  }
  col.gs <- c('pink', 'skyblue')
  names(col.gs) <- c('No', 'Yes')
  col.limmaPb <- c('pink', 'skyblue')
  names(col.limmaPb) <- c('nonDiff', 'Diff')
  annotation_colors[['gs']] <- col.gs
  annotation_colors[['limmaPb']] <- col.limmaPb
  
  col.signalType <- brewer.pal(8, 'Set3')[1:3]
  names(col.signalType) <- c('trend only', 'mean only', 'both')
  annotation_colors[['signalType']] <- col.signalType
  
  #### plot
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  if (break.0){
    cpl <- c(cpl[1:40], cpl[60:100])
  } 
  
  plist <- list()
  
  if (!is.na(sep)){
    rownames(expr.scale) <-sub(sep, '', rownames(expr.scale))
    rownames(rowann) <- sub(sep, ':.*', rownames(rowann))
    rownames(oridata) <- sub(sep, '', rownames(oridata))
  }
  
  rowann$meanDiff <- ifelse(rowMeans(oridata[sub(':.*','',rownames(rowann)),]) > 0,'Positive','Negative')
  
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
  colann.fit1 <-data.frame(pseudotime = rep(1:ncol(fit[[1]]), length(fit)),
                           group = gsub(sub('_.*', '_', names(fit)[1]),'',sub(';.*', '', colnames(fit.scale))), 
                           expression = 'ModelFitted',
                           stringsAsFactors = F)
  colann.fit2 <-data.frame(pseudotime = seq(1, ncol(FitDiff.scale)),
                           group = 'NA', 
                           expression = 'ModeledGroupDiff',
                           stringsAsFactors = F)
  colann.fit <- rbind(colann.fit1, colann.fit2)
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime)))
  names(col.pseudotime) = unique(colann.fit$pseudotime)
  annotation_colors$pseudotime <- col.pseudotime
  
  fit.scale <- cbind(fit.scale, FitDiff.scale)
  rownames(colann.fit) = colnames(fit.scale)
  
  if (!is.na(sep)){
    rownames(fit.scale) <-sub(sep, '', rownames(fit.scale))  
  }
  
  p2 <- pheatmap(
    fit.scale[, rownames(colann.fit)[colann.fit[,'expression'] != 'ModeledGroupDiff']],
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
    border_color = NA, silent = TRUE)
  plist[[3]] <- p2[[4]] 
  
  p3 <- pheatmap(
    fit.scale[, rownames(colann.fit)[colann.fit[,'expression'] == 'ModeledGroupDiff']],
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
    border_color = NA, silent = TRUE)
  
  p4data <- fit.scale[, rownames(colann.fit)[colann.fit[,'expression'] == 'ModeledGroupDiff']]
  p4 <- pheatmap(
    (p4data-rowMeans(p4data))/apply(p4data,1,sd),
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
    border_color = NA, silent = TRUE)
  
  
  plist[[5]] <- p3[[4]] 
  plist[[7]] <- p4[[4]]
  plist[[2]] <- plist[[4]] <- plist[[6]] <- ggplot(data=NULL) + geom_blank() + theme_void()
  # png(paste0(pdir, 'tmptest.png'),width = 11000,height = 10000,res = 300)
  grid.arrange(grobs = plist,layout_matrix=matrix(c(1,1,1,1,2,3,3,3,3,4,5,5,5,5,6,7,7,7,7),nrow=1))
  # dev.off()
}  
