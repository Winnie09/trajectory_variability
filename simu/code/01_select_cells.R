rm(list=ls())
setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
clu = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/cluster/resolution0.1.rds')
clu = as.numeric(as.character(clu))
u = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate/umap/umap.rds')
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
saveRDS(ord, '/Users/wenpinhou/Dropbox/trajectory_variability/simu/data/setting/select_cells_order.rds')

pat <- as.numeric(sub('BM','',sub(':.*','',sub('_.*','',n))))
names(pat) <- n
cluord <- order(sapply(1:6,function(i) mean(ord[names(tmpclu[tmpclu==i])])))
n <- names(tmpclu)
tmpclu <- match(tmpclu,cluord)
names(tmpclu) <- n
ord.bak = ord
source('/Users/wenpinhou/Dropbox/trajectory_variability/function/01_function.R')
f_vec <- sapply(seq(1,1e4), function(myseed){
  print(myseed)
  set.seed(myseed)
  v = runif(6)
  prop1 = v/sum(v)
  v = runif(6)
  prop2 = v/sum(v)
  prop <- prop1
  selectcell1 <- sapply(1:4,function(i) {
    tab = table(tmpclu[names(ord[pat==i])])
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
  prop <- prop2
  selectcell2 <- sapply(5:8,function(i) {
    tab = table(tmpclu[names(ord[pat==i])])
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
  # saveRDS(c(selectcell1, selectcell2), '/Users/wenpinhou/Dropbox/trajectory_variability/simu/data/setting/select_cells.rds')
  selectcell = c(selectcell1, selectcell2)
  ##########
  # rm(list=ls())
  # ord = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/simu/data/setting/select_cells_order.rds')
  # selectcell <- readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/simu/data/setting/select_cells.rds')
  ord = ord.bak[unlist(selectcell)]
  order = data.frame(Pseudotime = ord, Cell = names(ord), Patient = gsub(':.*','', names(ord)))
  g1 = names(selectcell)[1:4]
  g2 = names(selectcell)[5:8]
  f <- f_statistics_from_pseudotime(order, g1, g2)
})
  
saveRDS(f_vec, '/Users/wenpinhou/Dropbox/trajectory_variability/simu/data/setting/f_statistics.rds')

dist <- sapply(seq(1,1e4), function(myseed){
  set.seed(myseed)
  v = runif(6)
  prop1 = v/sum(v)
  v = runif(6)
  prop2 = v/sum(v)
  sum((prop1 - prop2)^2)
})
saveRDS(dist, '/Users/wenpinhou/Dropbox/trajectory_variability/simu/data/setting/prop_euclidean_dist.rds')


