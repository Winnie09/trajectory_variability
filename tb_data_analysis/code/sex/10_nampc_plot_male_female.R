library(here)
setwd(here())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
a = read.csv('tb/data/raw/mRNA_NAM_PCs.csv')
b = read.csv('tb/data/raw/CCA_NAM_PCs.csv')
c = readRDS('tb/data/sex/cellanno.rds')
design = readRDS('tb/data/sex/design.rds')

nampc1 = b[,4]
names(nampc1) = b[,1]
nampc2 = a[,4]
names(nampc2) = a[,3]

cell = sub('.*:', '', c[,1])
rownames(c) = cell

sex = design[,'sex']
names(sex) = sub('_.*', '', rownames(design))
sex = sex[sub('_.*', '', c[cell,2])]
sex = ifelse(sex == 1, 'female', 'male')

df = data.frame(nampc1 = nampc1[cell], nampc2 = nampc2[cell], sex = sex, stringsAsFactors = FALSE)

library(ggplot2)
library(scattermore)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
#source('resource/ggplot_theme.R')
theme_set(.new_theme)

pdir = 'tb/plot/sex/'
pdf(paste0(pdir, 'nampc_plot_male_female.pdf'), width = 2.7, height = 1.8)
ggplot(data = df, aes(x = nampc1, y = nampc2, color = sex)) + geom_scattermore(size = 0.001, stroke = 0) +
  scale_colour_manual(values = c('orange', 'royalblue')) +
  xlab('NAM PC1') + ylab('NAM PC2')+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()


