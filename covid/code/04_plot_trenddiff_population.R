Res = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/trenddiff.rds')
names(Res)
res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange, pvalue = Res$pvalue)
res <- res[res[,1] < 0.05, ]
rownames(res) <- sub(':.*', '', rownames(res))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- Res$pseudotime
expr = Res$expression
pred = Res$predict.values
rownames(expr) <- sub(':.*', '', rownames(expr))
Res$expression = expr
rownames(pred) = sub(':.*', '', rownames(pred))
Res$predict.values = pred
names(Res$parameter) = sub(':.*', '', names(Res$parameter))
names(Res$knotnum) = sub(':.*', '', names(Res$knotnum))
dg = rownames(res[res[,1]<0.05, ])

colnames(Res$design) <- c('intercept', 'condition')

library(splines)
fit <- sapply(dg, function(g){
  print(g)
  fit <- get_population_fit(Res, 'condition', g = g)
  vn <- sapply(1:length(fit), function(i){
    paste0(names(fit)[i], ';', rownames(fit[[i]]))
  })
  v <- as.vector(do.call(cbind, fit))
  names(v) <- vn
  v
})

fit = t(fit)

pd = reshape2::melt(fit)
colnames(pd) <- c('gene','cell','expr')
pd$covariate = sub(';.*', '', pd$cell)
pd$covariate[pd$covariate == 'condition_1'] = 'healthy'
pd$covariate[pd$covariate == 'condition_0'] = 'COVID19'
pd$t <- as.numeric(pseudotime[sub('.*;', '', pd[,2])])

library(ggplot2)
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/'
pdf(paste0(pdir, 'trenddiff_gene_population_curve.pdf'), width = 14, height = 8)
ggplot(data = pd, aes(x = t, y = expr, group = covariate,color = pd$covariate)) +
  geom_smooth() +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  xlab('Pseudotime') + ylab('Expression') +
  labs(color = NULL) +
  facet_wrap(~gene, scales = 'free') 
dev.off()
