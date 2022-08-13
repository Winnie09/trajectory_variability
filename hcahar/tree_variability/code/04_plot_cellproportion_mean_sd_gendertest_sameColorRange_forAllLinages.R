rm(list=ls())
library(RColorBrewer)
library(here)
library(ggplot2)
library(pheatmap)

## seurat 
res.ser <- read.csv(paste0('/home/whou10/scratch16/whou10/trajectory_variability/tree_variability/auto_pc_auto_nclu_module_3traj/result/sample.cellcomp.mean.csv'), row.names = 1)
rownames(res.ser) <- c('HSC->myeloid','HSC->erythroid','HSC->lymphocyte')
res.ser = res.ser[, c(2:5, 1, 6:8)]

## scvi
dir.scvi <- paste0(here(), '/hcaScvi/tree_variability/10latent_7clu/')
res.scvi <- read.csv(paste0(dir.scvi,'./res/sample.cellcomp.mean.csv'), row.names = 1)
rownames(res.scvi) <- c('HSC->myeloid','HSC->erythroid','HSC->lymphocyte')
res.scvi = res.scvi[, c(2:5, 1, 6:8)]

## harmony
dir.har <- paste0(here(), '/hcahar/tree_variability/')
res.har <- read.csv(paste0(dir.har,'./res/sample.cellcomp.mean.csv'), row.names = 1)
rownames(res.har) <- c('HSC->myeloid','HSC->erythroid','HSC->lymphocyte')
res.har = res.har[, c(2:5, 1, 6:8)]


min = range(cbind(res.ser, res.scvi, res.har))[1]
max = range(cbind(res.ser, res.scvi, res.har))[2]
breaksList = seq(min, max, by = 0.01)

gender <- data.frame(Sex=ifelse(1:8 %in% 1:4,'Male','Female'))
rownames(gender) <- colnames(res)

pdf(paste0(dir.har, './plot/cellproportion_mean_sd_gendertest_har.pdf'), width = 4, height= 4)
pheatmap(
  t(res.har),
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = gender,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
  breaks = breaksList
)
dev.off()

pdf(paste0(dir.har, './plot/cellproportion_mean_sd_gendertest_scvi.pdf'), width = 4, height= 4)
pheatmap(
  t(res.scvi),
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = gender,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
  breaks = breaksList
)
dev.off()


pdf(paste0(dir.har, './plot/cellproportion_mean_sd_gendertest_ser.pdf'), width = 4, height= 4)
pheatmap(
  t(res.ser),
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = gender,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
  breaks = breaksList
)
dev.off()


