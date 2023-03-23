setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/hscMepEry')
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/simu/data/hscMepEry')
pr = readRDS('./pca/pr.rds')
res <- UMAP(pr, save.umap=T, result.dir=getwd())

