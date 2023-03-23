library(here)
setwd(here())
source('function/01_function.R')
df = data.frame(scvi = c('erythroid', 'myeloid', 'lymphoid'), 
                seu = c('erythroid', 'monocyte', 'lymph'))
ovlist <- ov.pmlist <- list()
fdrmat = matrix(0, nrow = 2, ncol = 3)
rownames(fdrmat) = c('seurat', 'scVI')
colnames(fdrmat) = df[,1]

for (j in 1:3){
  print(j)
  print(df[j, ])
  scvi.res = readRDS(paste0('hcaScvi/testtime/res/', df[j, 1], '.rds'))
  seu.res = readRDS(paste0('hca/real/testtime/result/EM_pm/', df[j, 2], '/testtime_res.rds'))

  seu.s = seu.res[[1]]
  seu.s = as.data.frame(seu.s)
  seu.s = seu.s[order(seu.s[1], -seu.s[3]), ]
  
  scvi.s = scvi.res[[1]]
  scvi.s = as.data.frame(scvi.s)
  scvi.s = scvi.s[order(scvi.s[1], -scvi.s[3]), ]  
  
  int <- intersect(rownames(scvi.s), rownames(seu.s))
  str(int)
  scvi.s <- scvi.s[rownames(scvi.s) %in% int, ]
  seu.s <- seu.s[rownames(seu.s) %in% int, ]
  
  ## use Seurat significant genes as gs
  seu.s <- seu.s[seu.s[,1] < 0.05, ]
  str(seu.s)
  ## calculate overlap
  ovlist[[j]] = sapply(1:floor(nrow(scvi.s)/100), function(i){
    o = intersect(rownames(seu.s), rownames(scvi.s)[1:(i*100)])
    length(o)/nrow(seu.s)
  })  
  ovlist[[j]][floor(nrow(scvi.s)/100) + 1] <- length(intersect(rownames(seu.s), rownames(scvi.s)))/nrow(seu.s)
  
  set.seed(1234)
  scvi.sp = sample(rownames(scvi.s))
  
  ov.pmlist[[j]] = sapply(1:floor(length(scvi.sp)/100), function(i){
    o = intersect(rownames(seu.s), scvi.sp[1:(i*100)])
    length(o)/nrow(seu.s)
  })  
  ov.pmlist[[j]][floor(length(scvi.sp)/100) + 1] <- length(intersect(rownames(seu.s), scvi.sp))/nrow(seu.s)
  
  fdrmat[1 ,j] = nrow(seu.s)
  fdrmat[2 ,j] = sum(scvi.s[,1] < 0.05)
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
v = sapply(v, function(i) ifelse(i> nrow(scvi.s), nrow(scvi.s), i))
rownames(ov.pm) <- rownames(ov) <- paste0('top', v)
colnames(ov.pm) <- colnames(ov) <- df[,1]
  
  
saveRDS(ov, 'hcaScvi/testtime/compare/overlap_with_seurat_sig_res.rds')
saveRDS(ov.pm, 'hcaScvi/testtime/compare/overlap_with_seurat_sig_res_pm.rds')
saveRDS(fdrmat, 'hcaScvi/testtime/compare/fdrmat_seurat_sig_and_harmony.rds')

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
pdf('hcaScvi/testtime/plot/overlap_with_seurat_sig.pdf', height = 1.5, width = 2.4)
ggplot(data = pd3, aes(x = num, y = prop, color = lineage)) + 
  geom_line(aes(linetype = type), size = 0.3, alpha = 0.9) + scale_color_brewer(palette = 'Dark2') +
  theme(axis.text.x = element_text(angle=45, vjust = 0.5)) + 
  xlab('Number of top genes') + 
  ylab('Overlap proportion') 
dev.off()



