library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
lin = c('erythroid', 'monocyte', 'lymph')

ov <- ov.pm <- list()
fdrmat = matrix(0, nrow = 2, ncol = length(lin))
rownames(fdrmat) = c('seurat', 'harmony')
colnames(fdrmat) = lin

for (j in lin){
  print(j)
  har.res = readRDS(paste0('hcahar/real/testvar/result/', j, '/EM_pm.rds'))
  seu.res = readRDS(paste0('hca/real/testvar/result/EM_pm/', j, '/gender_res.rds'))
  har.s = har.res[[1]]
  har.s = as.data.frame(har.s)
  har.s = har.s[order(har.s[1], -har.s[3]), ]  
  
  seu.s = seu.res[[1]]
  seu.s = as.data.frame(seu.s)
  seu.s = seu.s[order(seu.s[1], -seu.s[3]), ]
  
  ov[[j]] = sapply(1:50, function(i){
    o = intersect(rownames(seu.s)[1:(i*100)], rownames(har.s)[1:(i*100)])
    length(o)/(i*100)
  })  
  
  set.seed(12345)
  seu.sp = sample(rownames(seu.s))
  har.sp = sample(rownames(har.s))
  
  ov.pm[[j]] = sapply(1:50, function(i){
    o = intersect(seu.sp[1:(i*100)], har.sp[1:(i*100)])
    length(o)/(i*100)
  })  
  
  fdrmat[1 ,j] = sum(seu.s[,1] < 0.05)
  fdrmat[2 ,j] = sum(har.s[,1] < 0.05)
}



  
ov = do.call(cbind, ov)  
ov.pm = do.call(cbind, ov.pm)
colnames(ov.pm) <- colnames(ov) <- df[,1]  
rownames(ov.pm) <- rownames(ov) <- paste0('top', (1:50)*100)
saveRDS(ov, 'hcahar/real/testvar/compare/overlap_with_seurat_res.rds')
saveRDS(ov.pm, 'hcahar/real/testvar/compare/overlap_with_seurat_res_pm.rds')
saveRDS(fdrmat, 'hcahar/real/testvar/compare/fdrmat_seurat_and_harmony.rds')

pd = reshape2::melt(ov)
pd[,1] = as.numeric(sapply(pd[,1], function(i) gsub('top', '', i))) 
colnames(pd) = c('num', 'lineage', 'prop')
pd$type = 'test results'

pd2 = reshape2::melt(ov.pm)
pd2[,1] = as.numeric(sapply(pd2[,1], function(i) gsub('top', '', i))) 
colnames(pd2) = c('num', 'lineage', 'prop')
pd2$type = 'permuted'

pd3 = rbind(pd, pd2)
pd3[,4] = factor(pd3[,4], levels = rev(sort(unique(pd3[,4]))))
 
library(ggplot2)
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('hcahar/real/testvar/plot/overlap_with_seurat.pdf', height = 1.5, width = 2.4)
ggplot(data = pd3, aes(x = num, y = prop, color = lineage)) + 
  geom_point(aes(shape = type), size = 0.2) + 
  geom_line(aes(linetype = type), size = 0.1, alpha = 0.5) + scale_color_brewer(palette = 'Dark2') +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5)) + 
  xlab('number of top genes') + 
  ylab('overlap proportion') 
dev.off()


