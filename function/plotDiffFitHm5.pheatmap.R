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
# colann = NULL
# rowann = NULL
# annotation_colors = NULL
# subsampleCell = TRUE
# numSubsampleCell=1e3
# sep = ':.*'
# break.0 = TRUE

plotDiffFitHm5.pheatmap <- function(testobj, showRowName = FALSE, cellWidthTotal = 250, cellHeightTotal = 400, showCluster = FALSE, colann = NULL, rowann = NULL, annotation_colors = NULL, subsampleCell = TRUE, numSubsampleCell=1e3, sep = NA, break.0 = TRUE){
  ## cellHeightTotal: when showRowName = TRUE, cellHeightTotal is suggested to be ten times the number of genes (rows).
  ## showCluster: (no implemented yet). if TRUE, "cluster" should be a slot in testobj, and it will be label in the heatmap. If FALSE, no need to pass in "cluster". 
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(ggplot2)
  testvar = testobj$testvar
  fit <- testobj$populationFit
  if ('XDEType' %in% names(testobj)) {
    XDEType <- testobj$XDEType
  } else {
    XDEType <- getXDEType(testobj)
  }
  XDEType <- XDEType[rownames(testobj$covariateGroupDiff)]
  
  if (subsampleCell){
    id <- round(seq(1, ncol(fit[[1]]), length.out = numSubsampleCell))
    for (i in 1:length(fit)){
      fit[[i]] <- fit[[i]][, id]
    }
    if (sum(XDEType == 'meanSig',na.rm=T) > 0){
      meanid <- which(XDEType == 'meanSig')
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
    if (sum(XDEType == 'meanSig', na.rm = T) > 0){
      meanid <- which(XDEType == 'meanSig')
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
  
  max <- apply(abs(testobj$covariateGroupDiff), 1, max)
  alluniformdiff <- testobj$covariateGroupDiff/max
  
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
  
  
  fit.scale <- do.call(cbind, fit)
  fit.scale <- fit.scale[names(testobj$cluster), ]
  fit.scale <- scalematrix(fit.scale)
  colnames(fit.scale) <- paste0(rep(names(fit), each = ncol(fit.scale)/length(fit)), ';cell', seq(1, ncol(fit.scale)))
  
  changepoint <- sapply(names(clu), function(i){
    ap <- which(FitDiff.sd[i, -ncol(FitDiff.sd)] * FitDiff.sd[i, -1] < 0)
    ap[which.min(abs(ap-ncol(FitDiff.sd)/2))]
  })
  res <- data.frame(clu = clu, 
                    cor = sapply(names(clu), function(i) cor(FitDiff.scale[i, seq(1, ncol(FitDiff.scale))], seq(1, ncol(FitDiff.scale)))),
                    # changepoint = sapply(names(clu), function(i) which.min(abs(FitDiff.sd[i, seq(1, ncol(FitDiff.scale))]))),
                    changepoint = changepoint,
                    XDEType = XDEType[names(clu)])
  
  
  #res <- res[order(res$XDEType, res$clu, res$changepoint), ]
  res1 <- res[res$XDEType=='trendSig',]
  res2 <- res[res$XDEType=='bothSig',]
  res3 <- res[res$XDEType=='other',]
  res4 <- res[res$XDEType=='meanSig',]
  
  o1 <- rownames(res1)[order(res1$clu,res1$cor>0,res1$changepoint)]
  
  pn <- rowMeans(alluniformdiff[rownames(res2),,drop=FALSE])
  o2 <- rownames(res2)[order(pn > 0,res2$clu,res2$cor>0,res2$changepoint)]
  
  o3 <- rownames(res3)[order(res3$clu,res3$cor>0,res3$changepoint)]
  #o1 <- rownames(res1)[order(match(res1$clu,names(sort(tapply(res1$cor,res1$clu,mean)))),res1$cor > 0,res1$changepoint)]
  # o2 <- rownames(res2)[order(match(res2$clu,names(sort(tapply(res2$cor,res2$clu,mean)))),res2$cor > 0,res2$changepoint)]
  # o3 <- rownames(res3)[order(match(res3$clu,names(sort(tapply(res3$cor,res3$clu,mean)))),res3$cor > 0,res3$changepoint)]
  o4 <- rownames(res4)[order(res4$clu)]
  #res <- res[c(o1,o2,o4,o3), ] ## put others in the last
  res <- res[c(o1,o2,o4), ] ## put others in the last
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
    expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(testobj$design)[testobj$design[, testvar] == sub('.*_','', i)], 1]]
  })
  expr.scale <- do.call(cbind, tmp)
  expr.scale <- scalematrix(expr.scale)
  expr.scale <- expr.scale[rownames(fit.scale), ]  
  
  ### annotate rows and columns
  if (is.null(colann)){
    colann <- data.frame(
      # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
      pseudotime = testobj$pseudotime[colnames(expr.scale)],
      group = as.character(testobj$design[cellanno[match(colnames(expr.scale), cellanno[,1]),2], testvar]),
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
      XDEType = as.character(XDEType[names(clu)]),
      stringsAsFactors = F)
    rownames(rowann) = names(clu)
  }
  rowann <- rowann[rownames(fit.scale), ,drop=F]
  rowann[,'XDEType'] <- factor(as.character(rowann[,'XDEType']), levels = c('trendSig','meanSig','bothSig','nonDDG','other'))
  
  if (length(unique(clu)) < 8){
    col.clu = brewer.pal(8, 'Set1')[1:length(unique(clu))]
  } else {
    col.clu = colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(clu)))
    
  }
  names(col.clu) = unique(clu)
  
  if (is.null(colann)| is.null(annotation_colors)){
    
    col.meanDiff = c('blue','red')
    names(col.meanDiff) <- c('Positive','Negative')
    #col.XDEType = brewer.pal(8, 'Set3')[1:5]
    #names(col.XDEType) = c('trendSig','meanSig','bothSig','nonDDG','other')
    
    col.XDEType = brewer.pal(8, 'Set3')[1:3]
    names(col.XDEType) = c('trendSig','meanSig','bothSig')
    
    annotation_colors = list(
      pseudotime = col.pseudotime,
      group = col.group,
      expression = col.expression,
      cluster = col.clu,
      XDEType = col.XDEType)
      #meanDiff = col.meanDiff)
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
  
  #rowann$meanDiff <- ifelse(rowMeans(oridata[sub(':.*','',rownames(rowann)),]) > 0,'Positive','Negative')
  
  p1data <- expr.scale
  p1data[p1data > quantile(as.vector(p1data), 0.95,na.rm=T)] <- quantile(as.vector(p1data), 0.95,na.rm=T)
  p1data[p1data < quantile(as.vector(p1data), 0.05,na.rm=T)] <- quantile(as.vector(p1data), 0.05,na.rm=T)
  p1 <- pheatmap(
    p1data,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    #color=cpl,
    color = colorRampPalette(c('blue3','skyblue','white','pink','red3'))(50),
    breaks=seq(-max(abs(p1data)),max(abs(p1data)),length.out=50),
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
  
  p2data <- fit.scale[, rownames(colann.fit)[colann.fit[,'expression'] != 'ModeledGroupDiff']]
  p2data[p2data > quantile(as.vector(p2data), 0.99,na.rm=T)] <- quantile(as.vector(p2data), 0.99,na.rm=T)
  p2data[p2data < quantile(as.vector(p2data), 0.01,na.rm=T)] <- quantile(as.vector(p2data), 0.01,na.rm=T)
  p2 <- pheatmap(
    p2data,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    #color = cpl,
    color = colorRampPalette(c('blue3','skyblue','white','pink','red3'))(50),
    breaks=seq(-max(abs(p2data)),max(abs(p2data)),length.out=50),
    annotation_col = colann.fit,
    annotation_row = rowann,
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal*1.23 / ncol(fit.scale),
    cellheight = cellHeightTotal / nrow(fit.scale),
    border_color = NA, silent = TRUE)
  plist[[3]] <- p2[[4]] 
  
  p3data <- fit.scale[, rownames(colann.fit)[colann.fit[,'expression'] == 'ModeledGroupDiff']]
  
  p3data[p3data > quantile(as.vector(p3data), 0.99,na.rm=T)] <- quantile(as.vector(p3data), 0.99,na.rm=T)
  p3data[p3data < quantile(as.vector(p3data), 0.01,na.rm=T)] <- quantile(as.vector(p3data), 0.01,na.rm=T)
  p3 <- pheatmap(
    p3data,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = colorRampPalette(c('blue3','skyblue','white','pink','red3'))(50),
    breaks=seq(-max(abs(p3data)),max(abs(p3data)),length.out=50),
    annotation_col = colann.fit,
    annotation_row = rowann,
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal*1.23 / ncol(fit.scale),
    cellheight = cellHeightTotal / nrow(fit.scale),
    border_color = NA, silent = TRUE)
  
  p4data <- fit.scale[, rownames(colann.fit)[colann.fit[,'expression'] == 'ModeledGroupDiff']]
  p4data <- (p4data-rowMeans(p4data))/apply(p4data,1,sd)
  
  p4data[p4data > quantile(as.vector(p4data), 0.99,na.rm=T)] <- quantile(as.vector(p4data), 0.99,na.rm=T)
  p4data[p4data < quantile(as.vector(p4data), 0.01,na.rm=T)] <- quantile(as.vector(p4data), 0.01,na.rm=T)
  
  p4 <- pheatmap(
    p4data,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = colorRampPalette(c('blue3','skyblue','white','pink','red3'))(50),
    breaks=seq(-max(abs(p4data),na.rm=T),max(abs(p4data),na.rm=T),length.out=50),
    annotation_col = colann.fit,
    annotation_row = rowann,
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal*1.23 / ncol(fit.scale),
    cellheight = cellHeightTotal / nrow(fit.scale),
    border_color = NA, silent = TRUE)
  
  # p5data <- fit.scale[, rownames(colann.fit)[colann.fit[,'expression'] == 'ModeledGroupDiff']]
  rownames(alluniformdiff) <- sub(':.*','',rownames(alluniformdiff))
  p5data <- rowMeans(alluniformdiff[rownames(fit.scale),,drop=FALSE]) %*% matrix(1,nrow=1,ncol=ncol(p4data))
  rownames(p5data) <- rownames(p4data)
  
  p5 <- pheatmap(
    p5data,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = colorRampPalette(c('blue3','skyblue','white','pink','red3'))(50),
    breaks=seq(-max(abs(p5data)),max(abs(p5data)),length.out=50),
    # annotation_col = colann.fit, ### add 20210616
    annotation_row = rowann,
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal*1.23 / ncol(fit.scale),
    cellheight = cellHeightTotal / nrow(fit.scale),
    border_color = NA, silent = TRUE)
  
  plist[[5]] <- p3[[4]] 
  plist[[7]] <- p4[[4]]
  plist[[9]] <- p5[[4]]
  plist[[2]] <- plist[[4]] <- plist[[6]] <- plist[[8]] <- ggplot(data=NULL) + geom_blank() + theme_void()
  # png(paste0(pdir, 'DiffFitHm5.png'),width = 8000,height = 3000,res = 300)
  grid.arrange(grobs = plist,layout_matrix=matrix(c(1,1,1,1,1,2,3,3,3,3,3,4,5,5,5,5,6,7,7,7,7,8,9,9,9,9),nrow=1))
  # dev.off()
}  


