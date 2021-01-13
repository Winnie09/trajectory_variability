library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

type = 'temra_testtime'
rdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type)
Res <- readRDS(here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type, 'res.rds'))
clu <- readRDS(paste0(rdir, '/cluster.rds'))

for (i in 1:max(clu)){
  print(paste0('cluster', i))
  res.go <- myGO(sub(':.*', '', names(clu[clu==i])), sub(':.*', '',rownames(Res$expr.ori)))
  res.go <- res.go[res.go$FDR<0.05, ]
  res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
  print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
  write.csv(res.go, paste0(rdir, '/GO_cluster', i, '.csv'))
}


