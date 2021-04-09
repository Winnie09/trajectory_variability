rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/addMultiSignalUsingExpr/result/nullsimu/'
# g = 'ABHD14B'
# Res = readRDS(paste0(rdir, 'EM_SelectKnots/testres.rds'))
# res = data.frame(Res$fdr, Res$foldchange, Res$pvalue)
# print(res[grep(g, rownames(res)),])
# res = res[order(res[,1], -res[,2]), ]
# res = res[res[,1]<0.05, ]
# write.csv(res, paste0(rdir, 'EM_SelectKnots/differential_genes.csv'))
# 
# num.fp <- nrow(res)
# names(num.fp) <- 'our_centered'

Res = readRDS(paste0(rdir, 'tradeSeq/testres.rds'))
v <- sapply(Res, function(res){
  res = res[order(res[,3], -res[,1]), ]
  # print(res[grep(g, rownames(res)),])
  res = res[res[,3]<0.05, ]
  nrow(res)
})
names(v) = paste0('tradeSeq_', names(v))
num.fp <- v


res = readRDS(paste0(rdir, 'tscan/testres.rds'))
# res = res[order(res[,3], -res[,2]), ]
res = res[res[,3]<0.05, ]
v = nrow(res)
names(v) = 'tscan'
num.fp = c(num.fp, v)

res = readRDS(paste0(rdir, 'monocle2/testres.rds'))
# res = res[order(res[,3], -res[,2]), ]
res = res[res[,3]<0.05, ]
v = nrow(res)
names(v) = 'monocle2'
num.fp = c(num.fp, v)


res = readRDS(paste0(rdir, 'monocle3/testres.rds'))
res = res[res[,5]<0.05, ]
# res = res[order(res[,5], -res[,4]), ]
v = nrow(res)
names(v) = 'monocle3'
num.fp = c(num.fp, v)

Res = readRDS(paste0(rdir, 'EM_pm/testres.rds'))
res = Res$statistics
# print(res[grep(g, rownames(res)),])
res = res[res[,1]<0.05, ]
v <- nrow(res)
names(v) <- 'Lamian'
num.fp = c(num.fp, v)

Res = readRDS(paste0(rdir, 'EM_chisq/testres.rds'))
res = Res$statistics
# print(res[grep(g, rownames(res)),])
res = res[res[,1]<0.05, ]
v <- nrow(res)
names(v) <- 'EM_chisq'
num.fp = c(num.fp, v)

num.fp <- sort(num.fp)
saveRDS(num.fp, paste0(rdir, 'perf/number_FP.rds'))



pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/addMultiSignalUsingExpr/plot/nullsimu/'
library(ggplot2)
pd = data.frame(method = names(num.fp), num.FP = num.fp, stringsAsFactors = F)
pd <- pd[pd$method != 'EM_chisq', ]
pd$method =factor(pd$method, levels = names(num.fp))
pdf(paste0(pdir, 'perf.pdf'),width=3.5,height=2)
ggplot() + geom_bar(data = pd, aes(x = method, y = num.FP, fill = method), stat = 'identity') +
  theme_classic() +
  coord_flip() +
  theme(legend.position = 'none') + 
  scale_fill_brewer(palette = 'Dark2') +
  theme(axis.text.y = element_text(color = c('red', rep('black', nrow(pd)-1)))) +
  ylab('number of false positives') + 
  xlab('method')
dev.off()

## plot
source('function/01_function.R')
plotGene(Res.EM, 'ABHD14B:ENSG00000114779')
plotGene(Res.EM, 'ABHD14B:ENSG00000114779', original.expr = T,facet.sample = T)


gg = "H1FX:ENSG00000184897"
pdf(paste0(pdir, sub(':.*', '', g), '_notcentered_point,pdf'), width = 3, height = 3)
plotGene(Res.center, gg, plot.point = T, point.size = 0.1, point.alpha = 0.5, original.expr = T, line.size = 0.5)
dev.off()
pdf(paste0(pdir, sub(':.*', '', g), '_notcentered_curve,pdf'), width = 4, height = 2.5)
plotGene(Res.center, gg,  original.expr = T, line.size = 0.5)
dev.off()

pdf(paste0(pdir, sub(':.*', '', g), '_center_point.pdf'), width = 4, height = 2.5)
plotGene(Res.center, gg, plot.point = T, point.size = 0.1, point.alpha = 0.5, original.expr = F, line.size = 0.5)
dev.off()
pdf(paste0(pdir, sub(':.*', '', g), '_center_curve.pdf'), width = 4, height = 2.5)
plotGene(Res.center, gg,  original.expr = F, line.size = 0.5)
dev.off()

