rm(list=ls())
library(RColorBrewer)
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/')
library(ggplot2)
library(pheatmap)
res <- read.csv(paste0('./auto_pc_auto_nclu_module_3traj/result/sample.cellcomp.mean.csv'), row.names = 1)
rownames(res) <- c('HSC->myeloid','HSC->erythroid','HSC->lymphocyte')
res = res[, c(2:5, 1, 6:8)]
sex <- data.frame(Sex=ifelse(1:8 %in% 1:4,'Male','Female'))
rownames(sex) <- colnames(res)
gender <- data.frame(Sex=ifelse(1:8 %in% 1:4,'Male','Female'))
rownames(gender) <- colnames(res)

############## branch proportion test
# test <- apply(res,1,function(i) t.test(i[2:5],i[-c(2:5)])$p.value)

# library(nnet)
# rs = reshape2::melt(as.matrix(res))
# colnames(rs) <- c('branch', 'sample', 'prop')
# # rs[,3] <- cut(rs[,3], breaks = seq(0, max(rs[,3])+0.1 ,0.1),include.lowest = TRUE)
# rs$sex <- sex[match(rs$sample, rownames(sex)),1]
# 
# rs$sex <- relevel(rs$sex, ref = "Male")
# test <- multinom(sex ~ prop + branch, data = rs)
# z <- summary(test)$coefficients/summary(test)$standard.errors
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
# p
library(nnet)
a <- readRDS('./auto_pc_auto_nclu_module_3traj/result/infer_tree_structure_res.rds')
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




##############

# df <- data.frame(mean=rowMeans(res),sd=apply(res,1,sd),pvalue=test, name=factor(rownames(res), levels = rownames(res)))

df <- data.frame(mean=rowMeans(res),sd=apply(res,1,sd), name=factor(rownames(res), levels = rownames(res)))
rownames(df) <- rownames(res)

pd = data.frame(t(res), Gender = gender, stringsAsFactors = F)
pd.point = reshape2::melt(t(res))
colnames(pd.point) = c('sample','name','prop')

write.csv(df, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/sourcedata/2D_bar.csv')

write.csv(pd, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/sourcedata/2D_hm.csv')


pdf('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/cellproportion_mean_sd_gendertest.pdf', width = 4, height= 4)
pheatmap(t(res),cluster_rows = F,cluster_cols = F,annotation_row = gender)
dev.off()

pdf('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/cellproportion_mean_sd_barplot.pdf', width = 2.6, height= 1.5)
ggplot()+
  geom_bar(data = df, aes(x = name, y = mean), stat="identity", fill=brewer.pal(11,'RdYlBu')[3], alpha=0.2) +
   geom_errorbar(data = df, aes(x=name, ymin=mean-sd, ymax=mean+sd), width=0.4, colour=brewer.pal(11,'RdYlBu')[10], alpha=0.9, size=1) +
  geom_jitter(data = pd.point, aes(x = name, y = prop), alpha = 0.3, size = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab('Sample proportion')  + xlab('')
dev.off()



