library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
source('function/01_function.R')
lin = c('erythroid', 'monocyte', 'lymph')
ovlist <- ov.pmlist <- list()
fdrmat = matrix(0, nrow = 2, ncol = length(lin))
rownames(fdrmat) = c('seurat', 'harmony')
colnames(fdrmat) = lin
for (j in lin){
  print(j)
  har.res = readRDS(paste0('hcahar/real/testvar/result/', j, '/EM_pm.rds'))
  seu.res = readRDS(paste0('hca/real/testvar/plot/EM_pm/', j, '/gender/numeric_res_with_clu.rds'))
  seu.s = seu.res[[1]]
  seu.s = as.data.frame(seu.s)
  seu.s = seu.s[order(seu.s[1], -seu.s[3]), ]
  
  har.s = har.res[[1]]
  har.s = as.data.frame(har.s)
  har.s = har.s[order(har.s[1], -har.s[3]), ]  
  
  int <- intersect(rownames(har.s), rownames(seu.s))
  str(int)
  har.s <- har.s[rownames(har.s) %in% int, ]
  seu.s <- seu.s[rownames(seu.s) %in% int, ]
  
  ## use Seurat significant genes as gs
  seu.s <- seu.s[rownames(seu.s) %in% names(seu.res[['cluster']]), ]
  str(seu.s)
  ## calculate overlap
  ovlist[[j]] = sapply(1:floor(nrow(har.s)/100), function(i){
    o = intersect(rownames(seu.s), rownames(har.s)[1:(i*100)])
    length(o)/nrow(seu.s)
  })  
  ovlist[[j]][floor(nrow(har.s)/100) + 1] <- length(intersect(rownames(seu.s), rownames(har.s)))/nrow(seu.s)
  
  set.seed(1234)
  har.sp = sample(rownames(har.s))
  
  ov.pmlist[[j]] = sapply(1:floor(length(har.sp)/100), function(i){
    o = intersect(rownames(seu.s), har.sp[1:(i*100)])
    length(o)/nrow(seu.s)
  })  
  ov.pmlist[[j]][floor(length(har.sp)/100) + 1] <- length(intersect(rownames(seu.s), har.sp))/nrow(seu.s)
  
  fdrmat[1 ,j] = nrow(seu.s)
  fdrmat[2 ,j] = sum(har.s[,1] < 0.05)
}
ovlist.bak = ovlist
n = max(sapply(ovlist, length))
ov <- sapply(1:length(ovlist), function(i){
  length(ovlist[[i]]) <- n
  ovlist[[i]]
})

ovpmlist.bak = ov.pmlist
n = max(sapply(ov.pmlist, length))
ov.pm <- sapply(1:length(ov.pmlist), function(i){
  length(ov.pmlist[[i]]) <- n
  ov.pmlist[[i]]
})

v <- (1:nrow(ov))*100
v = sapply(v, function(i) ifelse(i> nrow(har.s), nrow(har.s), i))
rownames(ov.pm) <- rownames(ov) <- paste0('top', v)
colnames(ov.pm) <- colnames(ov) <- lin
  
  
saveRDS(ov, 'hcahar/real/testvar/compare/overlap_with_seurat_sig_res.rds')
saveRDS(ov.pm, 'hcahar/real/testvar/compare/overlap_with_seurat_sig_res_pm.rds')
saveRDS(fdrmat, 'hcahar/real/testvar/compare/fdrmat_seurat_sig_and_harmony.rds')

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

source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('hcahar/real/testvar/plot/overlap_with_seurat_sig.pdf', height = 1.5, width = 2.4)
ggplot(data = pd3, aes(x = num, y = prop, color = lineage)) + 
  geom_line(aes(linetype = type), size = 0.5, alpha = 0.9) + scale_color_brewer(palette = 'Dark2') +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5)) + 
  xlab('number of top genes') + 
  ylab('overlap proportion') 
dev.off()



