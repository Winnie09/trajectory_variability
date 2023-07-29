library(here)
setwd(here())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
a = read.csv('tb_data_analysis/data/raw/mRNA_NAM_PCs.csv')
# b = read.csv('tb/data/raw/CCA_NAM_PCs.csv')
c = readRDS('tb_data_analysis/data/sex/cellanno.rds')
design = readRDS('tb_data_analysis/data/sex/design.rds')
# nampc1 = b[,4]
# names(nampc1) = b[,1]
nampc2 = a[,4]
names(nampc2) = a[,3]

cell = sub('.*:', '', c[,1])
rownames(c) = cell

sex = design[,'sex']
names(sex) = sub('_.*', '', rownames(design))
sex = sex[sub('_.*', '', c[cell,2])]
sex = ifelse(sex == 1, 'female', 'male')

u = readRDS('tb_data_analysis/data/sex/umap.rds')
rownames(u) = sub('.*:', '', rownames(u))
u = u[cell, ]

## pca 5, 6 seperate male and female
library(ggplot2)
library(scattermore)
library(viridis)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
# source('resource/ggplot_theme.R')
theme_set(.new_theme)

pdir = 'tb_data_analysis/plot/sex/'
pt = readRDS('tb_data_analysis/data/sex/ptpc2.rds')
names(pt) = sub('.*:', '', names(pt))

pdf(paste0(pdir, 'umap_plot_male_ptpc2.pdf'), width = 2.4, height = 1.6)
df = data.frame(nampc1 = u[,1], nampc2 = u[,2], sex = sex, pseudotime = pt[cell], stringsAsFactors = FALSE)
saveRDS(df, paste0(pdir, 'umap_plot_male_ptpc2_6A.rds'))
ggplot(data = df, aes(x = nampc1, y = nampc2, color = pseudotime, alpha = 0.5)) + geom_scattermore(size = 0.001, stroke = 0) +
  #scale_colour_manual(values = c('orange', 'royalblue')) +
  #guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_viridis(option = "A", direction = -1)+
  #labs(fill = 'pseudotime')+
  xlab('UMAP1') + ylab('UMAP2')
dev.off()



