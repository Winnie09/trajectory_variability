library(tradeSeq)
seed <- as.numeric(commandArgs(trailingOnly = T)[1])
partid <- as.numeric(commandArgs(trailingOnly = T)[2])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/condiments/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

pid <- 2

#expr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/expr.rds'))
#pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/pt.rds'))
#cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/cellanno.rds')
#design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/testdata/design.rds')

expr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds'))
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

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
expr <- expr[rowMeans(expr>0.1)>0.01, ]
pt <- pt[cellanno[,1]]
x <- design[match(cellanno[,2],rownames(design)),2]

x <- factor(x)
sce <- fitGAM(counts = expr, pseudotime = matrix(pt,ncol=1), cellWeights = matrix(1,nrow=ncol(expr),ncol=1),conditions=x, verbose = TRUE)
saveRDS(sce,file=paste0(rdir, seed,'_',partid,'sce_res.rds'))
cond_genes <- conditionTest(sce)
saveRDS(cond_genes,file=paste0(rdir, seed,'_',partid,'cond_gene_res.rds'))



