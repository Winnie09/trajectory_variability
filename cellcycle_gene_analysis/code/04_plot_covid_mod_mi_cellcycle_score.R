pd = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/res/covid_Mod_Mi_cellcycle_score.rds')
library(ggplot2)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)

library(RColorBrewer)
colv = brewer.pal(8, 'Dark2')[1:2]
names(colv) = c('Mi', 'Mod') 

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/plot/covid_Mod_Mi_cellcycle_score_s.pdf', width = 2, height = 1.5)
ggplot(data = pd) + geom_line(aes(x = pt.sample, y = cc.s.sample, color = severity, group = sample), size = 0.1, alpha = 0.8) + xlab('pseudotime') + ylab('cellcycle score (s genes)') + scale_color_manual(values = colv)
dev.off()

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/cellcycle/plot/covid_Mod_Mi_cellcycle_score_g2m.pdf', width = 2, height = 1.5)
ggplot(data = pd) + geom_line(aes(x = pt.sample, y = cc.g2m.sample, color = severity, group = sample), size = 0.1, alpha = 0.8) + xlab('pseudotime') + ylab('cellcycle score (G2-M genes)') + scale_color_manual(values = colv)
dev.off()




