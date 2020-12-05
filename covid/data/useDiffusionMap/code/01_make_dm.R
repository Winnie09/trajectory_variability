library(destiny)
library(here)
here()
seu <- readRDS(here('covid', 'data','GSE155673_pbmc','CD8integrate.rds'))
expr <- as.matrix(seu@assays$RNA@data)
str(expr)
expr <- expr[rowMeans(expr>0.1)>0.01, ] ## nothing changed
str(expr) ##
dm <- DiffusionMap(t(expr), )
str(dm)
saveRDS(dm, here('covid', 'data', 'useDiffusionMap', 'result','dm.obj.rds'))
dmap <- dm@@eigenvectors
saveRDS(dmap, here('covid','data','useDiffusionMap','result','dm.rds'))


pdf(here('covid','data','useDiffusionMap','plot', 'dm.pdf'), width = 6, height = 6)
plot(dm,1:2,
pch = 20) # pch for prettier points
dev.off()


