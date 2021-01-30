setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
order = readRDS('./hca/result/ery/order.rds')
mat = readRDS('./hca/data/HCA/proc/matrix/saver.rds')
mat = mat[,order$Cell]
source('./function/01_function.R')
order = data.frame(order, Patient = gsub('_.*','', order$Cell))
ap = as.character(unique(order$Patient))
age = as.numeric(sapply(ap, function(i) strsplit(i,':')[[1]][2]))
g1 = ap[order(age)[1:4]]
g2 = ap[order(age)[5:8]]
f_gene = f_statistics_from_gene(mat, order, g1, g2)
saveRDS(f_gene,'./hca/geneexpr/result/f_statistics_from_gene_age.rds')

mat <- mat[rowMeans(mat>0.01)>0.1, ]
a = f_statistics_from_gene_permute(mat, order, g1, g2, num.permute=1e4)
saveRDS(a,'./hca/geneexpr/result/f_statistics_from_lowExprGene_age_permute_new1e4.rds')
