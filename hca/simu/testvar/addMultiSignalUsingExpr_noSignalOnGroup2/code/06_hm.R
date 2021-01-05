library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))

for (dataType in seq(1, 4)) {
  print(dataType)
  ddir <-
    here('hca', 'data', 'simu', 'testvar', 'addMultiSignalUsingExpr')
  rdir <-
    here('hca',
         'simu',
         'testvar',
         'addMultiSignalUsingExpr',
         'result',
         'cluster',
         dataType)
  pdir <-
    here('hca',
         'simu',
         'testvar',
         'addMultiSignalUsingExpr',
         'plot',
         dataType)
  Res <-
    readRDS(
      here(
        'hca',
        'simu',
        'testvar',
        'addMultiSignalUsingExpr',
        'result',
        'EM_NOT_centered',
        paste0(dataType, '.rds')
      )
    )
  
  clu <- readRDS(paste0(rdir, '/cluster.rds'))
  fit <- readRDS(paste0(rdir, '/population_fit.rds'))
  selgene <- readRDS(paste0(ddir, '/selgene/selgene.rds'))
  
  allg <- names(sort(Res$fdr[Res$fdr < 0.05]))
  allg.p <- allg[1:min(length(allg), 25)]
  pdf(
    paste0(pdir, '/diffgene_population_fit.pdf'),
    width = 10,
    height = 10
  )
  plotGenePopulation(Res, allg.p, variable = 'group')
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
  
  
  pd.fit = t(sapply(fit, as.vector))
  fit.scale <-
    cbind(scalematrix(pd.fit[, 1:(ncol(pd.fit)/2)]),  scalematrix(pd.fit[, (ncol(pd.fit)/
                                                                                2 + 1):ncol(pd.fit)]))
  fit.scale <- fit.scale[names(clu), ]
  colnames(fit.scale) <- seq(1, ncol(fit.scale))
  
  res <- data.frame(clu = clu, 
                    cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(ncol(fit.scale)/2 + 1, ncol(fit.scale))], seq(ncol(fit.scale)/2 + 1, ncol(fit.scale)))))
    
  fit.scale <- fit.scale[rownames(res)[order(res$clu, res$cor)], ]
  
  
  cellanno <- Res$cellanno
  expr = Res$expr.ori
  expr <- expr[, names(Res$pseudotime)]
  expr.scale <-
    cbind(scalematrix(expr[allg, colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(Res$design[Res$design[, 2] ==
                                                                                                       0, ]), 1]]),
          scalematrix(expr[allg, colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(Res$design[Res$design[, 2] ==
                                                                                                       1, ]), 1]]))
  expr.scale <- expr.scale[rownames(fit.scale), ]
  ## plot ------------------------
  expr.scale[expr.scale > quantile(as.vector(expr.scale), 0.99)] <-
    quantile(as.vector(expr.scale), 0.99)
  expr.scale[expr.scale < quantile(as.vector(expr.scale), 0.01)] <-
    quantile(as.vector(expr.scale), 0.01)
  fit.scale[fit.scale > quantile(as.vector(fit.scale), 0.99)] <-
    quantile(as.vector(fit.scale), 0.99)
  fit.scale[fit.scale < quantile(as.vector(fit.scale), 0.01)] <-
    quantile(as.vector(fit.scale), 0.01)
  ### annotate rows and columns
  colann <-
    data.frame(
      sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
      pseudotime = Res$pseudotime[colnames(expr.scale)],
      addsignal = ifelse(Res$design[cellanno[match(colnames(expr.scale), cellanno[, 1]), 2], 2] ==
                           1, 'Yes', 'No'),
      stringsAsFactors = F
    )
  rownames(colann) = colnames(expr.scale)
  
  rowann = data.frame(
    cluster = as.character(clu),
    gs = sapply(names(clu), function(i)
      ifelse(i %in% selgene, 'Yes', 'No')),
    stringsAsFactors = F
  )
  rownames(rowann) = names(clu)
  #### define colors
  col.sample = colorRampPalette(rev(brewer.pal(n = 8, name = "Set1")))(length(unique(colann$sample)))
  names(col.sample) = unique(colann$sample)
  col.pseudotime = colorRampPalette(rev(brewer.pal(n = 9, name = "YlGnBu")))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  col.addsignal <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(colann$addsignal))]
  names(col.addsignal) <- sort(unique(colann$addsignal))
  col.clu = colorRampPalette(brewer.pal(8, 'Set1'))(max(clu))
  names(col.clu) = unique(clu)
  col.gs <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(rowann$gs))]
  names(col.gs) <- sort(unique(rowann$gs))
  #### save png
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  png(
    paste0(pdir, '/hm_kmeans_original_expression_scale.png'),
    width = 2100,
    height = 3200,
    res = 300
  )
  pheatmap(
    expr.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann,
    annotation_row = rowann,
    annotation_colors = list(
      sample = col.sample,
      pseudotime = col.pseudotime,
      addsignal = col.addsignal,
      cluster = col.clu,
      gs = col.gs
    ),
    cellwidth = 250 / ncol(expr.scale),
    cellheight = 450 / nrow(expr.scale),
    border_color = NA
  )
  dev.off()
  
  ## plot fitting values
  colann <-
    data.frame(addsignal = rep(c('Yes', 'No'), each = ncol(fit.scale) / 2),
               stringsAsFactors = F)
  rownames(colann) = colnames(fit.scale)
  
  col.addsignal <-
    brewer.pal(n = 8, name = "Pastel1")[1:length(unique(colann$addsignal))]
  names(col.addsignal) <- unique(colann$addsignal)
  
  png(
    paste0(pdir, '/hm_kmeans_population_fit_scale.png'),
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
      addsignal = col.addsignal,
      cluster = col.clu,
      gs = col.gs
    ),
    cellwidth = 250 / ncol(fit.scale),
    cellheight = 450 / nrow(fit.scale),
    border_color = NA
  )
  dev.off()
}
