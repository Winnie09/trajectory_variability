library(parallel)
library(VGAM)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/case/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/pt.rds')
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/case_design.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/expr.rds')

x <- design[match(cellanno[,2],rownames(design)),2]

pval <- mclapply(1:nrow(expr),function(i) {
  print(i)
  s <- expr[i,]
  full_model_fit <- VGAM::vglm(s~x+sm.ns(pt, df=3), epsilon=1e-1, family='uninormal')
  reduced_model_fit <- VGAM::vglm(s~x, epsilon=1e-1, family='uninormal')
  lrt <- VGAM::lrtest(full_model_fit,reduced_model_fit) 
  c(lrt@Body["Chisq"][2,],lrt@Body["Pr(>Chisq)"][2,])
},mc.cores=detectCores())
pval <- do.call(rbind,pval)
colnames(pval) <- c('stat','pval')
res <- data.frame(pval)
res$fdr=p.adjust(res$pval,method='fdr')
rownames(res) <- rownames(expr)
saveRDS(res, paste0(rdir,'monocle2_trajtest.rds'))  


