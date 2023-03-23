library(tradeSeq)
library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/monocyte/'
rdir <- paste0('hca/real/testvar/result/condiments/monocyte/gender/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
source('./function/01_function.R')

d = readRDS(paste0(ddir, 'input_expr.rds'))
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
pt = readRDS(paste0(ddir, 'input_pseudotime.rds'))

design = design[, 1:2]
design[,2] <- ifelse(design[,2] == 'male', 0, 1)
####

d <- d[,names(pt)]
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
x <- design[match(cellanno[,2],rownames(design)),2]
x <- factor(x)
sce <- fitGAM(counts = d, pseudotime = matrix(1:ncol(d),ncol=1), cellWeights = matrix(1,nrow=ncol(d),ncol=1),condition=x, verbose = TRUE)
saveRDS(sce,file=paste0(rdir, 'sce_res_forEff.rds'))
cond_genes <- conditionTest(sce)
saveRDS(cond_genes,file=paste0(rdir, 'cond_gene_res_forEff.rds'))


