library(here)
setwd(here())
source('function/01_function.R')
path <- 'erythroid'
ddir <- paste0('hca/real/build_from_tree_variability/result/', path, '/')
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/cellprop/data/')
dir.create(rdir, recursive = T)
cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
design = readRDS(paste0(ddir, 'input_design.rds'))
opt = readRDS(paste0(ddir, 'input_pseudotime.rds'))
design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
# library(ggplot2)
# ggplot(data.frame(pt=pt,samp=sub('_.*','',names(pt)),group=sub('_.*','',names(pt)) %in% rownames(design)[design[,2]==1]),aes(pt,color=group,group=samp)) + geom_density() + theme_classic()
for (seed in 1:1000) {
  print(seed)
  set.seed(seed)
  pt <- sample(opt)
  names(pt) <- names(opt)
  pt <- sort(pt)
  samp <- sub(':.*','',names(pt))
  names(samp) <- names(pt)
  for (prop in seq(0,0.5,0.05)) {
    spt <- pt
    exccell <- names(pt)[pt < median(pt)]
    spt <- spt[!names(spt) %in% sample(exccell,length(exccell)*prop)]
    saveRDS(spt,file=paste0(rdir, seed,'_',prop,'.rds'))
  }
}

