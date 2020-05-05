setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/cellcomposition/')
f = readRDS('./result/f_statistics_from_pseudotime_gender.rds')
f_p = data.frame(f_pm = readRDS('./result/f_statistics_from_pseudotime_gender_permute.rds'))
pval1 = sum(f_p>f)/nrow(f_p)

f2 = readRDS('./result/f_statistics_from_pseudotime_age.rds')
f2_p = data.frame(f_pm = readRDS('./result/f_statistics_from_pseudotime_age_permute.rds'))
pval2 = sum(f2_p>f2)/nrow(f2_p)

library(ggplot2)
library(gridExtra)

p1 <- ggplot() + geom_density(data=f_p, aes(x=f_pm), fill='pink', alpha=.3) + geom_vline(xintercept = f, col='red') + theme_classic() + ggtitle(paste0('Gender p=',pval1)) + 
  xlab('f statistics') + ylab('Density (1e4 Permutation)')


p2 <- ggplot() + geom_density(data=f2_p, aes(x=f_pm), fill='skyblue', alpha=.3) + geom_vline(xintercept = f2, col='darkblue') + theme_classic() + ggtitle(paste0('Age p=',pval2)) + 
  xlab('f statistics') + ylab('Density (1e4 Permutation)')


pdf('./plot/f_statistics.pdf',width=5,height=2.5)
grid.arrange(p1,p2,nrow=1)
dev.off()
