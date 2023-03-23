library(VGAM)
library(parallel)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/saverlog_pm.rds')
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/data/data/'
pt = readRDS(paste0(ddir, 'pseudotime_pm.rds'))
d <- d[,names(pt)]
design = cbind(1,matrix(c(1,1,0,0,1,1,0,0), nrow=8))
rownames(design) = paste0('BM',seq(1,8))
cellanno = data.frame(cell=colnames(d), sample = sub(':.*','', colnames(d)), stringsAsFactors = FALSE)
cellanno <- cellanno[match(colnames(d),cellanno[,1]),]
pt <- 1:ncol(d)
x <- design[match(cellanno[,2],rownames(design)),2]

pval <- mclapply(1:nrow(d),function(i) {
  s <- d[i,]
  full_model_fit <- VGAM::vglm(s~x+sm.ns(pt, df=3), epsilon=1e-1, family='uninormal')
  reduced_model_fit <- VGAM::vglm(s~x, epsilon=1e-1, family='uninormal')
  lrt <- VGAM::lrtest(full_model_fit,reduced_model_fit) 
  lrt@Body["Pr(>Chisq)"][2,]
},mc.cores=20)
pval <- unlist(pval)
res <- data.frame(pval=pval,fdr=p.adjust(pval,method='fdr'))
rownames(res) <- rownames(d)

saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/nullsimu_pm_window/result/monocle2_trajtest/res.rds')


