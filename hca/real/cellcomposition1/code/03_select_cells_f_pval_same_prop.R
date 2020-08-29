rm(list=ls())
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
clu = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/HCA/proc/cluster/resolution0.1.rds')
clu = as.numeric(as.character(clu))
u = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/HCA/proc/integrate/umap/umap.rds')
names(clu) <- row.names(u)

library(TSCAN)
selectid <- which(clu %in% c(7,14,8,3,0,4))
tmpclu <- clu[selectid]
n <- names(tmpclu)
tmpclu <- as.numeric(as.factor(tmpclu))
names(tmpclu) <- n
set.seed(12345)
ord <- rev(TSCANorder(exprmclust(t(u[selectid,]),cluster=tmpclu,reduce = F),orderonly = T))
n <- ord
ord <- 1:length(ord)
names(ord) <- n
saveRDS(ord, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/setting/select_cells_order.rds')

pat <- as.numeric(sub('BM','',sub(':.*','',sub('_.*','',n))))
names(pat) <- n
cluord <- order(sapply(1:6,function(i) mean(ord[names(tmpclu[tmpclu==i])])))
n <- names(tmpclu)
tmpclu <- match(tmpclu,cluord)
names(tmpclu) <- n
ord.bak = ord
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/setting/same_prop_tmp/')
f_pval <- sapply(seq(1,1e3), function(myseed){
  print(myseed)
  set.seed(myseed)
  v = runif(6)
  prop1 = v/sum(v)
  v = runif(6)
  prop2 = v/sum(v)
  prop <- prop1
  selectcell1 <- sapply(1:4,function(i) {
    
    tab = table(tmpclu[names(ord.bak[pat==i])])
    numCell = floor(min(tab/prop) * 0.8)
    set.seed(12345)
    for (dump in seq(1,1e4)){
      ct.alloc = table(sample(1:6, numCell, prob = prop, replace=T))  
      if (length(ct.alloc)==6) break   
    }
    unlist(sapply(1:6,function(j) {
      sample(names(ord[intersect(names(which(pat==i)), names(which(tmpclu==j)))]),min(length(names(ord[intersect(names(which(pat==i)), names(which(tmpclu==j)))])),ct.alloc[j]))
    }))
  })
  names(selectcell1) <- sapply(selectcell1, function(i) sub(':.*','',i[1]))
  # prop <- prop2
  selectcell2 <- sapply(5:8,function(i) {
    tab = table(tmpclu[names(ord.bak[pat==i])])
    numCell = floor(min(tab/prop) * 0.8)
    set.seed(12345)
    for (dump in seq(1,1e4)){
      ct.alloc = table(sample(1:6, numCell, prob = prop, replace=T))  
      if (length(ct.alloc)==6) break   
    }
    
    unlist(sapply(1:6,function(j) {
      sample(names(ord[intersect(names(which(pat==i)), names(which(tmpclu==j)))]),min(length(names(ord[intersect(names(which(pat==i)), names(which(tmpclu==j)))])),ct.alloc[j]))
    }))
  })
  names(selectcell2) <- sapply(selectcell2, function(i) sub(':.*','',i[1]))
  # saveRDS(c(selectcell1, selectcell2), '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/setting/select_cells.rds')
  selectcell = c(selectcell1, selectcell2)
  ord = ord.bak[unlist(selectcell)]
  order = data.frame(Pseudotime = ord, Cell = names(ord), Patient = gsub(':.*','', names(ord)))
  ## randomization 
  num.permute = 1e4
  
  f_permute <- sapply(seq(1, num.permute), function(myseed2){
    
    set.seed(myseed2)
    id = sample(1:nrow(order),nrow(order), replace=T) 
    names(selectcell) <- sample(names(selectcell))  
    ##
    g1 = names(selectcell)[1:4]
    g2 = names(selectcell)[5:8]
    f <- f_statistics_from_pseudotime(order[id,], g1, g2)
  })
  
  g1 = names(selectcell)[1:4]
  g2 = names(selectcell)[5:8]
  f <- f_statistics_from_pseudotime(order, g1, g2)
  result <- c(f,sum(f_permute >= f)/num.permute)
  saveRDS(result, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/setting/same_prop_tmp/',myseed,'.rds'))
  return(result)
})
rownames(f_pval) <- c('f','pvalue')
colnames(f_pval) <- paste0('seed_', seq(1,ncol(f_pval)))
saveRDS(f_pval, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/data/setting/f_pval_same_prop.rds')

fdr = p.adjust(f_pval[2,], method='fdr')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/simu/plot/setting_dist_f_pval_fdr_same_prop.pdf',width=7,height=3)
par(mfrow=c(1,2))
plot(f_pval[1,],f_pval[2,], main='', xlab='log10(f)', ylab='log10(p-value)', pch=20, cex=.5, log='xy')
plot(f_pval[1,],fdr, main='', xlab='log10(f)', ylab='log10(fdr)', pch=20, cex=.5, log='xy')
dev.off()
