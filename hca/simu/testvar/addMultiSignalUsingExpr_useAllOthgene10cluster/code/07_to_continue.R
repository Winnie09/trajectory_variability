# new a function
# input: testpt output including beta, phi
# input: covariables values user wants to know , if NULL then the unique values of the covatiates in the testpt data. if in the data only have age = 10, 20, 30, users can input 25 then we can output the pseudotime pattern of age == 25.
# phi * x * beta

signal = 1
method = 'EM'
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/perf/'
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testvar/addMultiSignalUsingExpr/selgene/selgene.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca')
rdir <- './simu/testvar/addMultiSignalUsingExpr/result/'
ddir <- './data/simu/testvar/addMultiSignalUsingExpr/'
r <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/EM/1.rds')
res = data.frame(fdr = r$fdr, foldchange = r$foldchange)
expr = r$expression
pseudotime = r$pseudotime
cellanno = r$cellanno

suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))

res <- res[order(res$fdr),]
tg <- rownames(res)[res$fdr < 0.05]
p <- sub(':.*','',colnames(expr))
expr <- expr[,order(p,pseudotime[colnames(expr)])]

set.seed(12345)
clu <- kmeans(expr[tg,],10,iter.max = 1000)$cluster
clu <- sort(clu)
saveRDS(clu, paste0(rdir, 'clu/geneclu.rds'))

colann <- data.frame(sample=sub(':.*','',colnames(expr)))
rownames(colann) <- colnames(expr)
colann$group <- as.factor(ifelse(colann[,1] %in% paste0('BM', c(3:4, 7:9)), 'Addsignal', 'Nosignal'))

rowann <- data.frame(cluster = clu[tg],
                     spikein = ifelse(tg %in% selgene, 'Yes', 'No'))
rownames(rowann) <- tg

library(pheatmap)
c1 <- brewer.pal(8,'Set1')
names(c1) <- paste0('BM',1:8)
c2 <- brewer.pal(8,'Dark2')[1:2]
names(c2) <- unique(colann$group)
r1 <- brewer.pal(10,'Set3')
names(r1) <- unique(rowann$cluster)
cpl <- colorRampPalette(rev(brewer.pal(n = 7, name =  "RdYlBu")))(100)
cpl <- c(rep(cpl[1],30),cpl,rep(cpl[length(cpl)],10))

# gaps_row=cumsum(rle(clu[1:30])$lengths)
# gaps_col = cumsum(rle(as.character(colann[,1]))$lengths)

png('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/clu/geneclu.30.10.png')
# png('geneclu.30.10.png')
pheatmap(expr[names(clu),],
              color=cpl, 
              cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
              annotation_col = colann,
              annotation_row = rowann[names(clu),,drop = F],
              annotation_colors = list(patient=c1, group = c2, cluster = r1))
dev.off()

###########################################

# data <- as.character(commandArgs(trailingOnly = T)[[1]][1])
# data = 'clusterType10_1'
# source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')

# plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/plot/', data, '/')
# rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/result/', data, '/')
# Res <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/result/', data, '/ptest_res.rds'))

Res = r ###

# res <- Res$res
# head(res)
# res <- res[res[,1] < 0.05 | res[,4] < 0.05, ]
res = data.frame(fdr = r$fdr, foldchange = r$foldchange)
res = res[order(res[,1], -res[,2]), ]
gene <- rownames(res[res[,1] < 0.05, ])
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


