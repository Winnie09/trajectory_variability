library(here)
setwd(here())
ddir <- 'genderdiff/data/'
library(data.table)

meta <- fread(paste0(ddir, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'), header = T, data.table = FALSE)
pheno <- fread(paste0(ddir, 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'), header = TRUE, data.table = FALSE)
rownames(pheno) <- pheno[,1]
meta <- meta[meta[,'SMTSD'] == 'Whole Blood', ]

m <- fread(paste0(ddir, 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'), data.table = FALSE)
m <- m[!duplicated(m[,2]), ]
rownames(m) <- m[,2]
m <- m[, c(-1,-2)]
m <- as.matrix(m)
m <- m[, colSums(m > 0) > 500]
m <- log2(m + 1)
m <- m[rowMeans(m > 0.1) > 0.01, ]
ap <- sapply(colnames(m), function(i) paste0(strsplit(i, '-')[[1]][1:2], collapse = '-'))
agg <- sapply(unique(ap), function(i){
  rowMeans(m[,ap == i, drop = FALSE])
})

age <- pheno[colnames(agg), 'AGE']
age.m <- sapply(age, function(i) mean(c(as.numeric(strsplit(i, '-')[[1]][1]), as.numeric(strsplit(i, '-')[[1]][2]))))
design = data.frame(intercept = 1, age = age.m)

library(limma)
res <- topTable(eBayes(lmFit(agg, design)),n=nrow(m),coef=2)
cor <- sapply(1:nrow(agg), function(i) cor(agg[i,], age.m))

names(cor) <- rownames(m)
head(rev(sort(abs(cor))))


for (path in c('monocyte', 'lymph', 'erythroid')){
  Res <- readRDS(paste0('hca/real/testvar/result/', path, '/age_fdr_res.rds'))
  Res <- Res[Res[,1] < 0.05, ]
  print(path)
  print(dim(Res))
  print(summary(abs(cor[sub(':.*', '', rownames(Res))])))
}
  
