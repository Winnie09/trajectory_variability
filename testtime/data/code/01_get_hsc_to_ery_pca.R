setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testtime/data/data/null/'
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
pr <- PCA(genebycell_mat = expr, save.pca = T, plot.statistics=T, plot.dir = rdir, result.dir = rdir, PC_for_columns = TRUE, findVariableGenes = TRUE, maxVariableGenes = 3000, numPC = NULL, smoothFittingMethod = 'loess')
UMAP(samplebyfeature_mat = pr$x[,1:23], save.umap = T, result.dir = rdir)
