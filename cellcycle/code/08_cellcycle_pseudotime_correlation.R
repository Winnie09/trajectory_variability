d  = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/res/tb_cellcycle_score.rds')
s = sub(':.*', '', rownames(d))
mat.tb = as.data.frame(t(sapply(unique(s), function(i){
  c(cor(d[s == i,1], d[s == i,3]),  cor(d[s == i,2], d[s == i,3]))
})))
colnames(mat.tb) = c('S', 'G2M')
library(reshape2)
pd = melt(mat.tb)
colnames(pd) = c('cellcycle', 'correlation')
pd$data = 'TB'


d  = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/res/covid_Mod_Mi_cellcycle_score.rds')

mat.covid = as.data.frame(t(sapply(unique(d[,4]), function(i){
  c(cor(d[d[,4] == i,1], d[d[,4] == i,3]),  cor(d[d[,4] == i,2], d[d[,4] == i,3]))
})))
colnames(mat.covid) = c('S', 'G2M')

pd2 = melt(mat.tb)
colnames(pd2) = c('cellcycle', 'correlation')
pd2$data = 'COVID19'

pd = rbind(pd, pd2)

library(ggplot2)
pd$data = as.factor(pd$data)
pd$grp = paste0(pd[,3], ';', pd[,1])
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/plot/cellcycle_pseudotime_correlation.pdf', width = 3.5, height = 2.4)
ggplot() + 
  geom_boxplot(data = pd, aes(x = grp, y = correlation, fill = cellcycle), outlier.shape = NA, position=position_dodge(1), alpha = 0.5)+
  scale_fill_brewer(palette = 'Set1')+
  geom_jitter(data = pd, aes(x = grp, y = correlation, color = cellcycle),  width = 0.25, size = 0.1, stoke = 0) + 
  ylab('correlation between cellcycle and density')
dev.off()  
  

