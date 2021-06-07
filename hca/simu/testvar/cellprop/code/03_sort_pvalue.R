setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- 'hca/simu/testvar/cellprop//result/'
dir.create(paste0(rdir, 'perf'), recursive = T)
reslist <- lapply(setdiff(list.files(rdir), 'perf'), function(path) {
  d = paste0(rdir, path,'/')
  df = lapply(list.files(d), function(f){
    readRDS(paste0(d, f))[[1]]
  })
  names(df) <- sub('.rds', '', sub('cell_proportion_test_', '',list.files(d)))
  do.call(rbind, df)
})
names(reslist) <- setdiff(list.files(rdir), 'perf')
saveRDS(reslist, paste0(rdir, 'perf/pvalue.rds'))


source('function/01_function.R')
rdir <- 'hca/simu/testvar/cellprop/result/'
for (path in setdiff(list.files(rdir),'perf')){
  d = paste0(rdir, path,'/')
  pdir <- paste0('hca/simu/testvar/cellprop/plot/', path,'/')
  dir.create(pdir, recursive = T)
  for (f in setdiff(list.files(d),'perf')){
    print(f)
    r = readRDS(paste0(d, f))
    pdf(paste0(pdir, sub('.rds','',f), '.pdf'), width = 4.5, height = 4)
    plotGene(r, 'prop',cellProp = T, variable = 'group', plot.point = T, facet.sample = T)
    dev.off()  
  } 
}
  
pdf('hca/simu/testvar/cellprop/plot/perf/pval.pdf', width = 6, height = 2)
par(mfrow = c(1,3))
for (i in 1:length(reslist)){
  plot(x = as.numeric(rownames(reslist[[i]])), y = reslist[[i]][,1], pch = 20, xlab = 'remove cell proportion', ylab='p.value', main = names(reslist)[i])
}
dev.off()
