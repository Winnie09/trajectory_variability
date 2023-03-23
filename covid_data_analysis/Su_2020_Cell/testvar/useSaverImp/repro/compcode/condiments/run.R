library(tradeSeq)
seed <- as.numeric(commandArgs(trailingOnly = T)[1])
partid <- as.numeric(commandArgs(trailingOnly = T)[2])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

rdir <- paste0('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/condiments/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
expr <- expr[, names(pt)]
expr <- expr[rowMeans(expr>0.1)>0.01, ]
expr <- expr[, cellanno[,1]]

set.seed(seed)
samp1 <- sample(rownames(design),nrow(design)/2)
samp2 <- setdiff(rownames(design),samp1)
cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
cellanno2 <- cellanno[cellanno[,2] %in% samp2,]
if (partid==1) {
  cellanno <- cellanno1
} else if (partid==2) {
  cellanno <- cellanno2
}

expr <- expr[,cellanno[,1]]
pt <- pt[cellanno[,1]]
x <- design[match(cellanno[,2],rownames(design)),2]

x <- factor(x)
sce <- fitGAM(counts = expr, pseudotime = matrix(pt,ncol=1), cellWeights = matrix(1,nrow=ncol(expr),ncol=1),conditions=x, verbose = TRUE)
saveRDS(sce,file=paste0(rdir, seed,'_',partid,'sce_res.rds'))
cond_genes <- conditionTest(sce)
saveRDS(cond_genes,file=paste0(rdir, seed,'_',partid,'cond_gene_res.rds'))


