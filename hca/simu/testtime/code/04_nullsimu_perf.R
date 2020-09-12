rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/nullsimu/'
g = 'ABHD14B'
Res = readRDS(paste0(rdir, 'EM_SelectKnots/testres.rds'))
res = data.frame(Res$fdr, Res$foldchange, Res$pvalue)
print(res[grep(g, rownames(res)),])
res = res[order(res[,1], -res[,2]), ]
res = res[res[,1]<0.05, ]
write.csv(res, paste0(rdir, 'EM_SelectKnots/differential_genes.csv'))

num.fp <- nrow(res)
names(num.fp) <- 'our_centered'


Res = readRDS(paste0(rdir, 'tradeSeq/testres.rds'))
v <- sapply(Res, function(res){
  res = res[order(res[,3], -res[,1]), ]
  print(res[grep(g, rownames(res)),])
  res = res[res[,3]<0.05, ]
  nrow(res)
})
names(v) = paste0('tradeSeq_', names(v))
num.fp <- c(num.fp, v)


res = readRDS(paste0(rdir, 'tscan/testres.rds'))
res = res[order(res[,3], -res[,2]), ]
res = res[res[,3]<0.05, ]
v = nrow(res)
names(v) = 'tscan'
num.fp = c(num.fp, v)

res = readRDS(paste0(rdir, 'monocle2/testres.rds'))
res = res[order(res[,3], -res[,2]), ]
res = res[res[,3]<0.05, ]
v = nrow(res)
names(v) = 'monocle2'
num.fp = c(num.fp, v)


res = readRDS(paste0(rdir, 'monocle3/testres.rds'))
res = res[res[,5]<0.05, ]
res = res[order(res[,5], -res[,4]), ]
v = nrow(res)
names(v) = 'monocle3'
num.fp = c(num.fp, v)

v = 0
names(v) = 'our_NOT_centered'
num.fp = c(num.fp, v)
num.fp = sort(num.fp)


library(ggplot2)
pd = data.frame(method = names(num.fp), num.FP = num.fp, stringsAsFactors = F)
pd$method =factor(pd$method, levels = names(num.fp))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/plot/nullsimu/perf.pdf',width=4,height=2)
ggplot() + geom_bar(data = pd, aes(x = method, y = num.FP, fill = method), stat = 'identity') +
  theme_classic() +
  coord_flip() +
  theme(legend.position = 'none')
dev.off()

## plot
plotGene(Res.EM, 'ABHD14B:ENSG00000114779')
plotGene(Res.EM, 'ABHD14B:ENSG00000114779', original.expr = T,facet.sample = T)
