library(Matrix)
source('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/function/01_function.R')

pid <- 2
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
expr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds'))
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

expr <- expr[rowMeans(expr>0.1)>0.01, ]

for (s in 1:10) {
  set.seed(s)
  samp1 <- sample(rownames(design),nrow(design)/2)
  samp2 <- setdiff(rownames(design),samp1)
  cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
  cellanno2 <- cellanno[cellanno[,2] %in% samp2,]
  
  res1 <- meandiff(expr = expr[,cellanno1[,1]], cellanno = cellanno1, design = design[samp1,])
  res2 <- meandiff(expr = expr[,cellanno2[,1]], cellanno = cellanno2, design = design[samp2,])
  
  saveRDS(list(res1,res2), paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/limma/',s,'.rds'))
}


