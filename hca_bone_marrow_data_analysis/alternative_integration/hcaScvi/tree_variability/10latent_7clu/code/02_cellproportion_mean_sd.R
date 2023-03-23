rm(list=ls())
library(RColorBrewer)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcaScvi/tree_variability/10latent_7clu/')
library(ggplot2)
library(pheatmap)
res <- read.csv(paste0('./res/sample.cellcomp.mean.csv'), row.names = 1)
rownames(res) <- c('HSC->myeloid','HSC->erythroid','HSC->lymphocyte')
res = res[, c(2:5, 1, 6:8)]
sex <- data.frame(Sex=ifelse(1:8 %in% 1:4,'Male','Female'))
rownames(sex) <- colnames(res)

library(nnet)
a <- readRDS('./res/infer_tree_structure_res.rds')
names(a)
o = a$order
names(o) <- paste0(gsub('.* ', 'c(', names(o)), ')')
k <- do.call(rbind, lapply(1:length(o), function(i){
  data.frame(cell = o[[i]], branch = names(o)[i])
}))

k$sample <- gsub(':.*', '', k$cell)
k$sex <- sex[match(k$sample, rownames(sex)), 1]
str(k)
test <- multinom(branch ~ sex, data = k)
z <- summary(test)$coefficients/summary(test)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

df <- data.frame(mean=rowMeans(res),sd=apply(res,1,sd), name=factor(rownames(res), levels = rownames(res)))
rownames(df) <- rownames(res)
gender <- data.frame(Sex=ifelse(1:8 %in% 1:4,'Male','Female'))
rownames(gender) <- colnames(res)

pdf('./plot/cellproportion_mean_sd_gendertest.pdf', width = 4, height= 4)
pheatmap(t(res),cluster_rows = F,cluster_cols = F,annotation_row = gender)
dev.off()

pdf('./plot/cellproportion_mean_sd_barplot.pdf', width = 2.6, height= 1.5)
ggplot(df)+
  geom_bar(aes(x = name, y = mean), stat="identity", fill=brewer.pal(11,'RdYlBu')[3], alpha=0.2) +
   geom_errorbar( aes(x=name, ymin=mean-sd, ymax=mean+sd), width=0.4, colour=brewer.pal(11,'RdYlBu')[10], alpha=0.9, size=1) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab('Sample proportion')  + xlab('')
dev.off()


