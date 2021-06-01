library(here)
setwd(here())
source('function/01_function.R')
for (path in c('erythroid', 'monocyte', 'lymph')){
  ddir <- paste0('hca/real/build_from_tree_variability/result/', path, '/')
  rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/data/',path,'/')
  dir.create(rdir, recursive = T)
  cellanno = readRDS(paste0(ddir, 'input_cellanno.rds'))
  design = readRDS(paste0(ddir, 'input_design.rds'))
  pt = readRDS(paste0(ddir, 'input_pseudotime.rds'))
  
  design = matrix( c(rep(1, 8),1,1,0,0,1,1,0,0), nrow=8)
  dimnames(design) = list(paste0('BM',seq(1,8)), c('intercept','group'))
  samp <- sub(':.*','',names(pt))
  names(samp) <- names(pt)
  set.seed(12345)
  for (prop in c(0.01,0.05, 0.25,seq(0.1, 0.9, 0.1))) {
    print(paste0(path, '_', prop))
    spt <- pt
    exccell <- intersect(names(samp)[which(samp %in% rownames(design)[design[,2]==1])],names(pt)[pt < median(pt)])
    spt <- spt[!names(spt) %in% sample(exccell,length(exccell)*prop)]
    saveRDS(spt,file=paste0(rdir, prop,'.rds'))
  }
}  


# ggplot(data.frame(pt=pt,samp=sub('_.*','',names(pt))),aes(pt,col=samp)) + geom_density() + facet_wrap(~samp)

# ggplot(data.frame(pt=spt,samp=sub('_.*','',names(spt))),aes(pt,col=samp)) + geom_density() + facet_wrap(~samp)

