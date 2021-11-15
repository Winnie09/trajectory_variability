library(VGAM)
library(parallel)
library(here)
setwd(here())
rdir <- paste0('covid/Su_2020_Cell/testvar/useSaverImp/result/monocle2_trajtest/Mod_Mi/')
dir.create(rdir, recursive = T, showWarnings = F)
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
d <- expr[rowMeans(expr>0.1)>0.01, ]

####

cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
pt <- seq(1,ncol(d))
x <- design[match(cellanno[,2],rownames(design)),2]

pval <- mclapply(1:nrow(d),function(i) {
print(i)
  s <- d[i,]
  full_model_fit <- VGAM::vglm(s~x+sm.ns(pt, df=3), epsilon=1e-1, family='uninormal')
  reduced_model_fit <- VGAM::vglm(s~x, epsilon=1e-1, family='uninormal')
  lrt <- VGAM::lrtest(full_model_fit,reduced_model_fit) 
  lrt@Body["Pr(>Chisq)"][2,]
},mc.cores=20)
pval <- unlist(pval)
res <- data.frame(pval=pval,fdr=p.adjust(pval,method='fdr'))
rownames(res) <- rownames(d)
saveRDS(res,file=paste0(rdir,'res.rds'))





