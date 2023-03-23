library(parallel)
library(VGAM)
seed <- as.numeric(commandArgs(trailingOnly = T)[1])
partid <- as.numeric(commandArgs(trailingOnly = T)[2])
print(seed)
print(partid)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')

rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/monocle2_trajtest/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

pid <- 2
expr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds'))
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

expr <- expr[rowMeans(expr>0.1)>0.01, ]

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

pval <- mclapply(1:nrow(expr),function(i) {
  print(i)
  s <- expr[i,]
  full_model_fit <- VGAM::vglm(s~x+sm.ns(pt, df=3), epsilon=1e-1, family='uninormal')
  reduced_model_fit <- VGAM::vglm(s~x, epsilon=1e-1, family='uninormal')
  lrt <- VGAM::lrtest(full_model_fit,reduced_model_fit) 
  c(lrt@Body["Chisq"][2,],lrt@Body["Pr(>Chisq)"][2,])
},mc.cores=10)
pval <- do.call(rbind,pval)
colnames(pval) <- c('stat','pval')
res <- data.frame(pval)
res$fdr=p.adjust(res$pval,method='fdr')
rownames(res) <- rownames(expr)
saveRDS(res, paste0(rdir, seed,'_',partid,'.rds'))  





