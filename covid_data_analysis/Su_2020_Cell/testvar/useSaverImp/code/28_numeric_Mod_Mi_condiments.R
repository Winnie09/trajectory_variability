library(tradeSeq)
library(here)
setwd(here())
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/condiments/Mod_Mi/')
dir.create(rdir, showWarnings = F, recursive = T)

d <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
pt <- readRDS('covid/Su_2020_Cell/data/tActivate_pseudotime.rds')
meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
cellanno <- readRDS('covid/Su_2020_Cell/data/cellanno.rds')
design <- readRDS('covid/Su_2020_Cell/data/design_numeric_Mod_Mi.rds')

rownames(cellanno) <- cellanno[,1]
cellanno <- cellanno[pt, ]
cellanno <- cellanno[cellanno[,2] %in% rownames(design), ]
pt <- pt[names(pt) %in% cellanno[,1]]
d <- d[, names(pt)]
d <- d[rowMeans(d>0.1)>0.01, ]

####
d <- d[,names(pt)]
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
x <- design[match(cellanno[,2],rownames(design)),2]
x <- factor(x)
sce <- fitGAM(counts = d, pseudotime = matrix(1:ncol(d),ncol=1), cellWeights = matrix(1,nrow=ncol(d),ncol=1),condition=x, verbose = TRUE)
saveRDS(sce,file=paste0(rdir, 'sce_res.rds'))
cond_genes <- conditionTest(sce)
saveRDS(cond_genes,file=paste0(rdir, 'cond_gene_res.rds'))



