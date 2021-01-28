rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# source('/Users/wenpinhou/Dropbox/resource/myfunc/01_function.R')
source('function/01_function.R')
for (dataType in seq(1, 4)) {
  print(dataType)
  ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr'
  rdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/cluster/',dataType)
  pdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/plot/', dataType)
  Res <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/EM_NOT_centered/',dataType,'.rds'))
  names(Res)
  # clu <- readRDS(paste0(rdir, '/cluster.rds'))
  fit <- readRDS(paste0(rdir, '/population_fit.rds'))
  str(fit)
  Res$populationFit <- fit
  fit.bak = fit
  # fit = fit[,1:(ncol(fit.bak)/2)]
  selgene <- readRDS(paste0(ddir, '/data/selgene/selgene.rds'))
  selgene1 <- readRDS(paste0(ddir, '/data/selgene/selgene1.rds'))
  selgene2 <- readRDS(paste0(ddir, '/data/selgene/selgene2.rds'))
  selgene3 <- readRDS(paste0(ddir, '/data/selgene/selgene3.rds'))
  
  allg <- names(sort(Res$fdr[Res$fdr < 0.05]))
  allg.p <- allg[1:min(length(allg), 25)]
  pdf(
    paste0(pdir, '/diffgene_population_fit.pdf'),
    width = 10,
    height = 10
  )
  plotGenePopulation(Res, allg.p, type = 'variable')
  dev.off()
  pdf(paste0(pdir, '/diffgene_by_sample.pdf'),
      width = 10,
      height = 10)
  plotGene(Res, allg.p)
  dev.off()
  pdf(paste0(pdir, '/diffgene_by_group.pdf'),
      width = 10,
      height = 10)
  plotGene(Res, allg.p, variable = 'group')
  dev.off()
  

  clu <- readRDS(paste0(rdir, '/cluster.rds'))
  str(clu)
  str(fit)
  # fit.scale <- lapply(1:length(fit), function(i){
  #   m <- scalematrix(fit[[i]])
  #   m <- m[names(clu), ]
  #   colnames(m) <- seq(1, ncol(m))
  #   m
  # })
  # names(fit.scale) <- names(fit)
  # str(fit.scale)
  # sd.scale <- lapply(1:length(fit), function(i){
  #   apply(fit[[i]], 1, sd)
  # })
  # str(sd.scale)
  # 
  # fit.scale <- do.call(cbind, fit.scale)
  
  fit.scale <- do.call(cbind, fit)
  fit.scale <- fit.scale[names(clu), ]
  fit.scale <- scalematrix(fit.scale)
  str(fit.scale)
  colnames(fit.scale) <- seq(1, ncol(fit.scale))
  
  res <- data.frame(clu = clu, 
                    cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(1, ncol(fit.scale)/2)], seq(1, ncol(fit.scale)/2))))
  res$signalType <- sapply(rownames(res), function(i) {
    if (i %in% selgene1){
      '1'
    } else if (i %in% selgene2){
      '2'
    } else {
      '3'
    }
  })
  fit.scale <- fit.scale[rownames(res)[order(res$clu, res$signalType, res$cor)], ]
  # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
  meanres <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/meandiff/', dataType, '.rds'))
  str(meanres)
  meanDiffTest <- ifelse(meanres[rownames(fit.scale), 5] < 0.05, 'Diff', 'nonDiff')
  names(meanDiffTest) <- rownames(fit.scale)
  ## ------------------------
  ## plot original expression 
  ## ------------------------
  cellanno <- Res$cellanno
  expr = Res$expr.ori
  expr <- expr[, names(Res$pseudotime)]
  # expr.scale <-
  #   cbind(scalematrix(expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(Res$design[Res$design[, 2] == 0, ]), 1]]),
  #         scalematrix(expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(Res$design[Res$design[, 2] == 1, ]), 1]]))
  expr.scale <-
    cbind(expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(Res$design[Res$design[, 2] == sub('.*_','',names(fit)[1]), ]), 1]],
          expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(Res$design[Res$design[, 2] == sub('.*_','',names(fit)[2]), ]), 1]])
  
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
  colann <- data.frame(sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
      pseudotime = Res$pseudotime[colnames(expr.scale)],
      addsignal = ifelse(Res$design[cellanno[match(colnames(expr.scale), cellanno[, 1]), 2], 2] == 1, 'Yes', 'No'),
      stringsAsFactors = F)
  rownames(colann) = colnames(expr.scale)
  
  rowann = data.frame(
    cluster = as.character(clu),
    signalType = sapply(names(clu), function(i){
      if (i %in% selgene1){
        'trenddiff only'
      } else if (i %in% selgene2){
        'meandiff only'
      } else {
        'trenddiff_meandiff'
      }
    }),
    meanDiffTest = meanDiffTest[names(clu)],
    gs = sapply(names(clu), function(i)
      ifelse(i %in% selgene, 'Yes', 'No')),
    stringsAsFactors = F
  )
  rownames(rowann) = names(clu)
  rowann <- rowann[rownames(fit.scale), ]
  #### define colors
  col.sample = colorRampPalette(rev(brewer.pal(n = 8, name = "Set1")))(length(unique(colann$sample)))
  names(col.sample) = unique(colann$sample)
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  col.addsignal <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(colann$addsignal))]
  names(col.addsignal) <- sort(unique(colann$addsignal))
  col.clu = colorRampPalette(brewer.pal(8, 'Set1'))(max(clu))
  names(col.clu) = unique(clu)
  col.signalType <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(rowann$signalType))]
  names(col.signalType) <- sort(unique(rowann$signalType))
  col.gs <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(rowann$gs))]
  names(col.gs) <- sort(unique(rowann$gs))
  col.meanDiffTest <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(rowann$meanDiffTest))]
  names(col.meanDiffTest) <- sort(unique(rowann$meanDiffTest))
  #### save png
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  png(
    paste0(pdir, '/hm_kmeans_original_expression_scale2.png'),
    width = 2100,
    height = 3200,
    res = 300
  )
  pheatmap(
    expr.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = F,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann,
    annotation_row = rowann,
    annotation_colors = list(
      sample = col.sample,
      pseudotime = col.pseudotime,
      addsignal = col.addsignal,
      cluster = col.clu,
      gs = col.gs,
      signalType = col.signalType,
      meanDiffTest = col.meanDiffTest
    ),
    cellwidth = 250 / ncol(expr.scale),
    cellheight = 450 / nrow(expr.scale),
    border_color = NA
  )
  dev.off()
  ## --------------------
  ## plot fitting values
  ## --------------------
  colann <-
    data.frame(pseudotime = rep(seq(1, ncol(fit.scale)/length(fit)), length(fit)),
               addsignal = rep(c('Yes', 'No'), each = ncol(fit.scale) / length(fit)),
               stringsAsFactors = F)
  rownames(colann) = colnames(fit.scale)
  
  col.addsignal <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(colann$addsignal))]
  names(col.addsignal) <- unique(colann$addsignal)
  
  png(
    paste0(pdir, '/hm_kmeans_population_fit_scale2.png'),
    width = 2100,
    height = 3200,
    res = 300
  )
  pheatmap(
    fit.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann,
    annotation_row = rowann,
    annotation_colors = list(
       pseudotime = col.pseudotime,
      addsignal = col.addsignal,
      cluster = col.clu,
      signalType = col.signalType,
      gs = col.gs,
      meanDiffTest = col.meanDiffTest
    ),
    cellwidth = 250 / ncol(fit.scale),
    cellheight = 450 / nrow(fit.scale),
    border_color = NA
  )
  dev.off()
}


str(res)
id <- rownames(res[res$clu == 1 & res$signalType == 3, ])
id <- id[10:30]
id
str(id)
g.tmp <- rownames(expr.scale[rownames(expr.scale) %in% id, ],)
str(g.tmp)
g = g.tmp[6]
g
plotGene(Res,g)
plotGenePopulation(Res, g, type = 'variable')

str(rowann)
id <- rownames(rowann)[which(rowann$gs == 'Yes' & rowann$meanDiffTest == 'nonDiff' & rowann$signalType == 'trenddiff only')]
str(id)

## plot trendonly signal & limman cannot identifieid genes
png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_genes.png'), width = 3000, height = 3000, res = 200)
plotGene(Res, id[1:100], plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
dev.off()
png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_genes_population.png'), width = 3000, height = 3000, res = 200)
plotGenePopulation(Res, id[1:100], type = 'variable', sep = ':.*')
dev.off()

## plot selected example genes
allg <- c('CTDP1', 'CFDP1', 'NOL6', 'TMEM198B')
allg <- sapply(allg, function(i){
  rownames(fit.scale)[grepl(i, rownames(fit.scale))]
})
png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_example.png'), width = 800, height = 550, res = 200)
plotGene(Res, allg, plot.point = T, variable = colnames(Res$design)[2], continuous = F, point.size = 0.1, point.alpha = 0.2, sep = ":.*")
dev.off()
png(paste0(pdir, '/meanDiffTest_notdiff_and_trenddiff_example_population.png'), width = 850, height = 500, res = 200)
plotGenePopulation(Res, allg, type = 'variable', sep = ':.*')
dev.off()

## plot cluster mean and difference
pdf(paste0(pdir, '/clusterMeanAndDiff.pdf'), width = 3, height = 6)
plotClusterMeanAndDiff(Res, clu)
dev.off()

