library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))

type = 'tex_testvar'
rdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type)
pdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','plot', type)
dir.create(pdir)
Res <- readRDS(here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type, 'res.rds'))
allg <- names(sort(Res$fdr[Res$fdr < 0.05]))
fit <- sapply(allg, function(g){
  tmp <- get_population_fit(Res, variable = 'type', gene = g)
}, simplify = FALSE)
saveRDS(fit, paste0(rdir, '/population_fit.rds')) ##########

if (length(allg) < 9){
  pdf(paste0(pdir, '/diffgene_population_fit.pdf'), width = 8, height = 8)
  print(plotGenePopulation(Res, allg, variable = 'type'))
  dev.off()
  pdf(paste0(pdir, '/diffgene_by_sample.pdf'), width = 6, height = 6)
  print(plotGene(Res, allg))
  dev.off()
  pdf(paste0(pdir, '/diffgene_by_group.pdf'), width = 6, height = 6)
  print(plotGene(Res, allg, variable = 'type'))
  dev.off()
}

pd.fit = t(sapply(fit, as.vector))
colnames(pd.fit) <- rep(seq(1,ncol(pd.fit)/2), 2)

cellanno <- Res$cellanno
expr = Res$expr.ori
expr <- expr[,names(Res$pseudotime)]
pd.expr <- cbind(scalematrix(expr[allg, colnames(expr) %in% cellanno[cellanno[,2] %in% rownames(Res$design[Res$design[,2]==0,]),1]]),
              scalematrix(expr[allg, colnames(expr) %in% cellanno[cellanno[,2] %in% rownames(Res$design[Res$design[,2]==1,]),1]]))

## plot
### annotate rows and columns
colann <-
data.frame(sample = cellanno[match(cellanno[, 1], colnames(pd.expr)), 2],
           pseudotime = Res$pseudotime[colnames(pd.expr)],
           stringsAsFactors = F)
rownames(colann) = colnames(pd.expr)

#### define colors
col.sample = colorRampPalette(rev(brewer.pal(n = 8, name = "Set1")))(length(unique(colann$sample)))
names(col.sample) = unique(colann$sample)
col.pseudotime = colorRampPalette(rev(brewer.pal(n = 9, name = "YlGnBu")))(length(unique(colann$pseudotime)))
names(col.pseudotime) = unique(colann$pseudotime)

library(ggplot2)

pd.expr[pd.expr > quantile(as.vector(pd.expr), 0.99)] <- quantile(as.vector(pd.expr), 0.99)
pd.expr[pd.expr < quantile(as.vector(pd.expr), 0.01)] <- quantile(as.vector(pd.expr), 0.01)
pd.fit[pd.fit > quantile(as.vector(pd.fit), 0.99)] <- quantile(as.vector(pd.fit), 0.99)
pd.fit[pd.fit < quantile(as.vector(pd.fit), 0.01)] <- quantile(as.vector(pd.fit), 0.01)


png(
  paste0(pdir, '/hm_kmeans_population_fit_scale.png'),
  width = 2100,
  height = 3200,
  res = 300
)
cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
pheatmap(
pd.expr,
cluster_rows = F,
cluster_cols = FALSE,
show_rownames = TRUE,
show_colnames = FALSE,
color = cpl,
annotation_col = colann,
annotation_colors = list(
  sample = col.sample
),
cellwidth = 250/ncol(pd.expr),
cellheight = 450/nrow(pd.expr),
border_color = NA
)
dev.off()


png(
  paste0(pdir, '/hm_kmeans_original_expression_scale.png'),
  width = 2100,
  height = 3200,
  res = 300
)
pheatmap(
pd.fit,
cluster_rows = F,
cluster_cols = FALSE,
show_rownames = TRUE,
show_colnames = FALSE,
color = cpl,
cellwidth = 250/ncol(pd.expr),
cellheight = 450/nrow(pd.expr),
border_color = NA
)
dev.off()




