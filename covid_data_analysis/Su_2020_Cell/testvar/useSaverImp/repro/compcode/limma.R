library(Matrix)
source('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/function/01_function.R')

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, cellanno[,1]]
expr <- expr[rowMeans(expr>0.1)>0.01, ]

for (s in 1:100) {
  set.seed(s)
  samp1 <- sample(rownames(design),nrow(design)/2)
  samp2 <- setdiff(rownames(design),samp1)
  cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
  cellanno2 <- cellanno[cellanno[,2] %in% samp2,]
  
  res1 <- meandiff(expr = expr[,cellanno1[,1]], cellanno = cellanno1, design = design[samp1,])
  res2 <- meandiff(expr = expr[,cellanno2[,1]], cellanno = cellanno2, design = design[samp2,])
  
  saveRDS(list(res1,res2), paste0('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/limma/Mod_Mi',s,'.rds'))
}

