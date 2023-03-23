# -----------------------------------------------
# clustering, population pattern, hm, GO analysis
# -----------------------------------------------
rm(list=ls())
method <- as.character(commandArgs(trailingOnly = T)[[1]])
print(method)
library(ggplot2)
library(reshape2)
library(scattermore)
library(RColorBrewer)
library(splines)
library(viridis)
library(here)
here()
source(here('function', '01_function.R'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
pdir <- here('covid', 'plot', method)
dir.create(pdir, recursive = T)
rdir <- here('covid', 'result', method, 'GO')
dir.create(rdir, recursive = T)
# clustering
Res <- readRDS(here('covid', 'result', method, paste0(method, '.rds')))
if (is.data.frame(Res)){
  if ('adj.P.Val' %in% colnames(Res)) diffgene <- rownames(Res)[Res[,'adj.P.Val'] < 0.05]
  if ('fdr' %in% colnames(Res)) diffgene <- rownames(Res)[Res[,'fdr'] < 0.05]
  if ('FDR' %in% colnames(Res)) diffgene <- rownames(Res)[Res[,'FDR'] < 0.05]
} else {
  if ('fdr' %in% names(Res)) diffgene <- names(Res$fdr[Res$fdr < 0.05])
}
  
if (length(diffgene) == 0)  {
  sink(paste0(rdir, '/readme.txt'))
  print('There is no differential genes!')
  sink()
} else {
  Res.more <- readRDS(here('covid', 'result', paste0(sub('_.*', '', method),'_trenddiff'), paste0(sub('_.*', '', method), '_trenddiff.rds')))
  pt <- Res.more$pseudotime
  pred <- Res.more$predict.values[diffgene, names(pt)]
  pred.scale <- (pred - rowMeans(pred))/apply(pred, 1, sd)
  set.seed(12345)
  clu <- kmeans(pred.scale, max(mykmeans(pred.scale)$cluster)+1)$cluster 
  names(clu) = rownames(pred)
  clu <- sort(clu)
  write.csv(data.frame(gene = names(clu), cluster = clu), paste0(rdir,'/', max(clu), 'cluster.csv'))
  
  pred <- pred[names(clu), ]
  pred.scale <- pred.scale[names(clu), ]
  
  agg <- t(sapply(unique(clu), function(i) colMeans(pred[clu == i, ,drop=F])))
  rownames(agg) <- paste0('cluster', unique(clu), '(',table(clu), ')')
  
  library(ggplot2)
  pd <- reshape2::melt(agg)
  pd$pseudotime = pt[pd[,2]]
  colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
  
  cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
  pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
  
  pdf(paste0(pdir,'/', max(clu), 'cluster_pattern_curve.pdf'), width = 5, height = 3.5)
  ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
          geom_smooth() +
          theme_classic() +
          scale_color_manual(values = colorRampPalette(brewer.pal(9,'Set1'))(max(clu))) +
          xlab('Pseudotime') +
          ylab('Averaged Fitted Centered Expression')
  dev.off()
  
  
  ##### hm
  # annotate rows and columns
  celltype <- readRDS(here('covid','data', 'allintegrate', 'celltype.rds'))
  cluster <- readRDS(here('covid','data', 'allintegrate', 'cluster.rds'))
  design <- readRDS(here('covid', 'data', 'tex', 'design.rds'))
  cellanno <- Res.more$cellanno
  colann <- data.frame(sample = cellanno[match(cellanno[,1], colnames(pred)),2],
                       stringsAsFactors = F)
  rownames(colann) = colnames(pred)
  rowann = data.frame(cluster = as.character(clu),stringsAsFactors = F)
  rownames(rowann) = names(clu)
  
  colann <- data.frame(colann, 
                       celltype = celltype[cluster[rownames(colann)]], 
                       gender = design[match(colann[,1], design[,6]),'gender'],
                       type = design[match(colann[,1], design[,6]),'type'],
                       stringsAsFactors = F)
  # define colors
  library(pheatmap)
  library(RColorBrewer)
  
  col.clu = colorRampPalette(brewer.pal(9,'Set1'))(max(clu))
  names(col.clu) = unique(clu)
  col.ct <- colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(colann$celltype)))
  names(col.ct) = unique(colann$celltype)
  col.ct = col.ct[!is.na(names(col.ct))]
  col.sample = colorRampPalette(rev(brewer.pal(n = 8, name = "Set1")))(length(unique(colann$sample)))
  names(col.sample) = unique(colann$sample)
  col.gender = colorRampPalette(brewer.pal(8,'Set2'))(length(unique(colann$gender)))
  names(col.gender) <- unique(colann$gender)
  col.type <- colorRampPalette(brewer.pal(8,'Accent'))(length(unique(colann$type)))
  names(col.type) <- unique(colann$type)
  
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  cpl = c(rep(cpl[1], 80), cpl, rep(cpl[length(cpl)], 120))
  png(paste0(pdir, 'hm_kmeans.png'),width=1050,height=1600, res=200)
  pheatmap(pred.scale[names(clu), ],
           cluster_rows = F, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE,
           color = cpl,
           annotation_row = rowann, annotation_col = colann,
           annotation_colors = list(cluster = col.clu, celltype = col.ct, sample = col.sample,
                                    gender = col.gender, type = col.type),
           cellwidth = 241.24/ncol(pred), cellheight = 467.35/nrow(pred))
  dev.off()
  
  # GO analysis  
  if (grepl('trenddiff', method)){
    fdr = Res$fdr
    fdr = fdr[order(fdr)]
    pdf(paste0(pdir, 'example_genes.pdf'),width=8,height=6.5)
    print(plotGene(Res, names(fdr[1:9])))
    dev.off()
  }
    
  for (i in 1:max(clu)){
    print(paste0('cluster', i))
    res.go <- myGO(names(clu[clu==i]), rownames(Res.more$expr.ori))
    res.go <- res.go[res.go$FDR<0.05, ]
    res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
    print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
    write.csv(res.go, paste0(rdir, '/cluster', i, '.csv'))
  }
}

  
