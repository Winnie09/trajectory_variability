smat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/hscMepEry/matrix/normcount.rds')
# ct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/hscMepEry/ct/ct.rds')
set.seed(12345)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/hscMepEry/pca/')
set.seed(12345)
pca = PCA(genebycell_mat=smat, save.pca = TRUE, plot.statistics=TRUE, plot.dir = getwd(), result.dir = getwd(), PC_for_columns = TRUE, findVariableGenes = TRUE, numPC = NULL)

