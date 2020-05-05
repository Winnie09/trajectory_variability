setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
order = readRDS('./hca/result/ery/order.rds')
mat = readRDS('./hca/data/HCA/proc/matrix/saver.rds')
mat = mat[,order$Cell]
source('./function/01_function.R')
order = data.frame(order, Patient = gsub('_.*','', order$Cell))
ap = as.character(unique(order$Patient))
g1 = ap[grepl('female', ap)]
g2 = ap[grepl(':male', ap)]
f_gene = f_statistics_from_gene(mat, order, g1, g2)
saveRDS(f_gene,'./hca/geneexpr/result/f_statistics_from_gene_gender.rds')

mat <- mat[rowMeans(mat>0.01)>0.1, ]
a = f_statistics_from_gene_permute(mat, order, g1, g2, num.permute=1e4)
saveRDS(a,'./hca/geneexpr/result/f_statistics_from_lowExprGene_gender_permute.rds')
