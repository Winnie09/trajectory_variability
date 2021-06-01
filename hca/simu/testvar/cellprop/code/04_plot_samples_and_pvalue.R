rdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/hca/simu/testvar/cellprop/result/'
for (path in list.files(rdir)){
  d = paste0(rdir, path,'/')
  pdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/hca/simu/testvar/cellprop/plot/', path,'/')
  dir.create(pdir, recursive = T)
  for (f in list.files(d)){
    r = readRDS(paste0(d, f))
    pdf(paste0(pdir, sub('.rds','',f), '.pdf'), width = 4.5, height = 4)
    plotGene(r, 'prop',cellProp = T, variable = 'group', plot.point = T, facet.sample = T)
    dev.off()  
  } 
}
  
pdf('/Users/wenpinhou/Dropbox/trajectory_variability/hca/simu/testvar/cellprop/plot/perf/pval.pdf', width = 6, height = 2)
par(mfrow = c(1,3))
for (i in 1:length(reslist)){
  plot(x = as.numeric(rownames(reslist[[i]])), y = reslist[[i]][,1], pch = 20, xlab = 'remove cell proportion', ylab='p.value', main = names(reslist)[i])
}
dev.off()

# r1 = readRDS(paste0(rdir, 'erythroid/cell_proportion_test_0.01.rds'))
# plotGene(r1, 'prop',cellProp = T, variable = 'group', plot.point = T, facet.sample = T)
# 
# r2 = readRDS(paste0(rdir, 'erythroid/cell_proportion_test_0.5.rds'))
# plotGene(r2, 'prop',cellProp = T, variable = 'group', plot.point = T, facet.sample = T)


