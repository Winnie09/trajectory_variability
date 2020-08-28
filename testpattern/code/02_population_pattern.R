# new a function
# input: testpt output including beta, phi
# input: covariables values user wants to know , if NULL then the unique values of the covatiates in the testpt data. if in the data only have age = 10, 20, 30, users can input 25 then we can output the pseudotime pattern of age == 25.
# phi * x * beta

data <- as.character(commandArgs(trailingOnly = T)[[1]][1])
# data = 'clusterType10_1'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/plot/', data, '/')
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/result/', data, '/')

Res <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/result/', data, '/ptest_res.rds'))
res <- Res$res
head(res)
res <- res[res[,1] < 0.05 | res[,4] < 0.05, ]
gene <- rownames(res[res[,4] < 0.05, ])
expr = Res$expression[gene, ]
knotnum = Res$knotnum[gene]
design = Res$design
cellanno = Res$cellanno
rownames(cellanno) = cellanno[,1]
pseudotime = Res$pseudotime

beta <- lapply(Res$trenddiff.parameter[gene], function(i){
  i$beta
})
names(beta) <- gene
predict.values <- Res$predict.values[gene,]

pseudotime = pseudotime[order(pseudotime)]
expr = expr[, names(pseudotime)]
predict.values = predict.values[, names(pseudotime)]
 
# plot(fit1~pseudotime, col = 'red', pch = 20, cex = .5,ylim=c(0,3))
# points(fit2~pseudotime, pch = 20, cex = .5)
library(splines)
fit <- sapply(gene, function(g){
  print(g)
  fit <- get_population_fit(Res, 'condition', g = g)
  vn <- sapply(1:length(fit), function(i){
    paste0(names(fit)[i], ';', rownames(fit[[i]]))
  })
  v <- as.vector(do.call(cbind, fit))
  names(v) <- vn
  v
})

clu <- mykmeans(fit,10)$cluster
agg <- aggregate(t(fit),list(clu),mean)
agg <- agg[,-1]
agg <- as.matrix(agg)
rownames(agg) <- paste0('cluster', seq(1, nrow(agg)))

saveRDS(fit, paste0(rdir, 'trenddiff_gene_popoulation_fit_10clu.rds'))
saveRDS(clu, paste0(rdir, 'trenddiff_gene_cluster_10clu.rds'))
saveRDS(agg, paste0(rdir, 'trenddiff_gene_agg_10clu.rds'))

library(reshape2)
pd <- melt(agg)
pd$covariate <- as.factor(sub(';.*', '', pd[,2]))
pd$x <- as.numeric(pseudotime[sub('.*;', '', pd[,2])])
pd$cell <- sub('.*;', '', pd$Var2)
library(ggplot2)

pdf(paste0(plotdir, 'population_level_trenddiff_gene_10clu.pdf'), width = 4, height = max(clu)*2)
ggplot() +
  geom_point(data = pd, aes(x = x, y = value, group = covariate,color = pd$covariate)) +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  xlab('Pseudotime') + ylab('Expression') +
  labs(color = NULL) +
  facet_wrap(~Var1, ncol = 1) 
dev.off()

library(pheatmap)
library(RColorBrewer)
hmpd <- predict.values[names(sort(clu)), names(pseudotime)]
hmpd <- hmpd[, order(as.character(pd[match(colnames(hmpd), pd$cell),'covariate']), pseudotime[colnames(hmpd)])]
anno_row <- data.frame(cluster = as.factor(clu[rownames(hmpd)]))
rownames(anno_row) = rownames(hmpd)
anno_col <- data.frame(covariate = pd[match(colnames(hmpd), pd$cell),'covariate'])
rownames(anno_col) <- colnames(hmpd)


png(paste0(plotdir, 'population_level_trenddiff_gene_hm_10clu.png'), width = 600, height = 600)
pheatmap(hmpd, cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_row = anno_row,
         annotaton_col = anno_col)
dev.off()


selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/data/data/selgene/selgene.rds')
selgene = sub(':.*','',selgene)
df <- data.frame(interceptdiff = ifelse(res[,1]<0.05, 'interceptdiff', 'nointerceptdiff'), trenddiff = ifelse(res[,4]<0.05, 'trenddiff','notrenddiff'), selgene = ifelse(rownames(res) %in% selgene, 'true', 'false'))
rownames(df) = rownames(res)

