library(VGAM)
library(parallel)
library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/lymph/'
rdir <- paste0('hca/real/testvar/result/monocle2_trajtest/lymph/gender/')
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


