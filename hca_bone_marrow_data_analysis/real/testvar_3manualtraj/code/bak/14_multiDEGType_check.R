library(here)
library(ggplot2)
setwd(here())
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
  
mat <- sapply(c('erythroid', 'lymph', 'monocyte'), function(path){
print(path)
res <- readRDS(paste0(rdir, path, '/gender_fdr_res.rds'))
res <- res[res[,7] < 0.05, ]
DDGType <- as.character(res[,10])
names(DDGType) <- rownames(res)
tab1 <- table(DDGType[rownames(res)])

res.lm <- readRDS(paste0(ddir, path,'/meandiff_gender_res.rds'))
res.lm <- res.lm[res.lm[5] < 0.05, ]
tab2 <- table(DDGType[rownames(res.lm)])
tab3 <- table(DDGType[setdiff(rownames(res), rownames(res.lm))])

mat <- matrix(0, nrow = 4, ncol=3)
dimnames(mat) <- list(names(tab1), c('ourmethod', 'limma', 'new'))
mat[names(tab1),1] <- tab1
mat[names(tab2),2] <- tab2
mat[names(tab3),3] <- tab3
mat
})

mat <- sapply(c('erythroid', 'lymph', 'monocyte'), function(path){
print(path)
res <- readRDS(paste0(rdir, path, '/gender_fdr_res.rds'))
res <- res[res[,7] < 0.05, ]
DDGType <- as.character(res[,10])
names(DDGType) <- rownames(res)
tab1 <- table(DDGType[rownames(res)])

res.lm <- readRDS(paste0(ddir, path,'/meandiff_gender_res.rds'))

res.lm <- res.lm[order(res.lm[,5], abs(res.lm[,1])), ]
res.lm <- res.lm[1:nrow(res), ]
tab2 <- table(DDGType[rownames(res.lm)])
tab3 <- table(DDGType[setdiff(rownames(res), rownames(res.lm))])

mat <- matrix(0, nrow = 4, ncol=3)
dimnames(mat) <- list(names(tab1), c('ourmethod', 'limma', 'new'))
mat[names(tab1),1] <- tab1
mat[names(tab2),2] <- tab2
mat[names(tab3),3] <- tab3
mat
})

par(mfrow = c(2,3))
for (path in c('erythroid', 'lymph', 'monocyte')){
print(path)
res <- readRDS(paste0(rdir, path, '/gender_fdr_res.rds'))
res.lm <- readRDS(paste0(ddir, path,'/meandiff_gender_res.rds'))
plot(res[,1], res.lm[rownames(res),5], xlab = 'ourmethod', ylab = 'limma', main = path)
}

