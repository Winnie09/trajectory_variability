library(tradeSeq)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/case/condiments/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/case_design.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/count.rds')

x <- design[match(cellanno[,2],rownames(design)),2]

x <- factor(x)
sce <- fitGAM(counts = expr, pseudotime = matrix(pt,ncol=1), cellWeights = matrix(1,nrow=ncol(expr),ncol=1),conditions=x, verbose = TRUE,BPPARAM = MulticoreParam(5))
saveRDS(sce,file=paste0(rdir,'sce_res.rds'))
cond_genes <- conditionTest(sce)
saveRDS(cond_genes,file=paste0(rdir, 'cond_gene_res.rds'))

