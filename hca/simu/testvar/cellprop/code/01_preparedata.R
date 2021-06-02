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
  g1g <- rownames(design)[design[,2]==1]
  g0g <- rownames(design)[design[,2]==0]
  
  window <- cut(pt,quantile(pt,seq(0,1,length.out=101)),include.lowest = T)
  names(window) <- names(pt)
  tab <- table(samp,window)
  prop <- tab/rowSums(tab)
  sn <- rowSums(tab)
  exclist <- unlist(sapply(unique(window),function(sw){
    g1swm <- median(prop[g1g,sw]) ## not median of median
    g0swm <- median(prop[g0g,sw])
    sid <- names(window)[window==sw]
    ssamp <- sub(':.*','',sid)
    if (g1swm < g0swm) { ## sub-sample cells from  each sample in the low group weighted by the sample's proportion in that group within this window (otherwise bugs when a sample does not have cells or enough cells in this window to sub-sample from)
      unlist(sapply(names(which(!is.na(table(ssamp)[g0g]))),function(si) sample(sid[ssamp==si],round(sn[si] * (g0swm - g1swm) * table(ssamp)[si]/length(ssamp) ))))
    } else {
      unlist(sapply(names(which(!is.na(table(ssamp)[g1g]))),function(si) {
        sample(sid[ssamp==si],round(sn[si] * (g1swm - g0swm) * table(ssamp)[si]/length(ssamp)))
        }))
    }
  }))
  
  pt <- pt[!names(pt) %in% exclist]
  saveRDS(pt,file=paste0(rdir, '0.rds'))
  for (prop in c(0.01,0.05, 0.25,seq(0.1, 0.9, 0.1))) {
    print(paste0(path, '_', prop))
    spt <- pt
    exccell <- intersect(names(samp)[which(samp %in% rownames(design)[design[,2]==1])],names(pt)[pt < median(pt)])
    spt <- spt[!names(spt) %in% sample(exccell,length(exccell)*prop)]
    saveRDS(spt,file=paste0(rdir, prop,'.rds'))
  }
}
