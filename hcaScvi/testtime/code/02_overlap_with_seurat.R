library(here)
setwd(here())
# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/testtime/res/')

df = data.frame(scvi = c('erythroid', 'myeloid', 'lymphoid'), 
                seu = c('erythroid', 'monocyte', 'lymph'))


ov <- ov.pm <- list()
fdrmat = matrix(0, nrow = 2, ncol = length(df[,1]))
rownames(fdrmat) = c('seurat', 'harmony')
colnames(fdrmat) = df[,1]

for (j in 1:3){
  print(df[j, ])
  scvi.res = readRDS(paste0('hcaScvi/testtime/res/', df[j, 1], '.rds'))
  seu.res = readRDS(paste0('hca/real/testtime/result/EM_pm/', df[j, 2], '/testtime_res.rds'))
  scvi.s = scvi.res[[1]]
  seu.s = seu.res[[1]]
  scvi.s = scvi.s[order(scvi.s[1], -scvi.s[3]), ]
  seu.s = seu.s[order(seu.s[1], -seu.s[3]), ]
  
  ov[[j]] = sapply(1:50, function(i){
    o = intersect(rownames(seu.s)[1:(i*100)], rownames(scvi.s)[1:(i*100)])
    length(o)/(i*100)
  })  
  
  set.seed(12345)
  seu.sp = sample(rownames(seu.s))
  scvi.sp = sample(rownames(scvi.s))
 
  ov.pm[[j]] = sapply(1:50, function(i){
    o = intersect(seu.sp[1:(i*100)], scvi.sp[1:(i*100)])
    length(o)/(i*100)
  })  
  
  fdrmat[1 ,j] = sum(seu.s[,1] < 0.05)
  fdrmat[2 ,j] = sum(scvi.s[,1] < 0.05)
}
ov = do.call(cbind, ov)  
ov.pm = do.call(cbind, ov.pm)
colnames(ov.pm) <- colnames(ov) <- df[,1]  
rownames(ov.pm) <- rownames(ov) <- paste0('top', (1:50)*100)
saveRDS(ov, 'hcaScvi/testtime/compare/overlap_with_seurat_res.res')
saveRDS(ov.pm, 'hcaScvi/testtime/compare/overlap_with_seurat_res_pm.res')
saveRDS(fdrmat, 'hcaScvi/testtime/compare/fdrmat_seurat_and_scvi.res')
write.csv(fdrmat, 'hcaScvi/testtime/plot/fdrmat_seurat_and_scvi.csv')


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

pdf('hcaScvi/testtime/plot/overlap_with_seurat.pdf', height = 1.4, width = 2.1)
ggplot(data = pd3, aes(x = num, y = prop, color = lineage)) + 
  geom_point(aes(shape = type) ,size = 0.2) + 
  geom_line(aes(linetype = type), size = 0.2, alpha = 0.5) + scale_color_brewer(palette = 'Dark2') +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5)) + 
  xlab('number of top genes') + 
  ylab('overlap proportion') 
dev.off()

