##  prepare data for differential test
### read in data, filtering
seu <- readRDS(here('covid', 'data','GSE155673_pbmc','integrated','CD8integrate.rds'))
expr <- as.matrix(seu@assays$RNA@data) ## librarysize-normalized  log2-transform
str(expr)
expr <- expr[rowMeans(expr>0.1)>0.01, ] ## nothing changed
str(expr) ##

ord1 <- readRDS(here('covid','data','GSE155673_pbmc','useDiffusionMap','tex','pseudotime_order.rds'))
ord2 <- readRDS(here('covid','data','GSE155673_pbmc','useDiffusionMap','temra','pseudotime_order.rds'))

### transform the cell names to be the same as SAVER imputed cell names
smeta <- read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/data/GSE155673_pbmc/meta/sampleMeta2.txt', stringsAsFactors = FALSE)
smeta[,1] <- sub('_', '', sub('_RNA', '', smeta[,1]))
rownames(smeta) <- smeta[,6]
colnames(smeta) <- c('samplename', 'type2','type3','age','gender','sample')
smeta$type <- sub('-.*', '', smeta[,2])
v1 <- sub(':.*', '', sub('GSE155673_', '', colnames(expr)))
v2 <- sub('.*:', '', sub('GSE155673_', '', colnames(expr)))
v1 <- smeta[match(v1,smeta[,1]),6]
colnames(expr) <- paste0(v1,':', v2)

v1 <- sub(':.*', '', sub('GSE155673_', '', ord1))
v2 <- sub('.*:', '', sub('GSE155673_', '', ord1))
v1 <- smeta[match(v1,smeta[,1]),6]
ord1 <- paste0(v1,':', v2)

v1 <- sub(':.*', '', sub('GSE155673_', '', ord2))
v2 <- sub('.*:', '', sub('GSE155673_', '', ord2))
v1 <- smeta[match(v1,smeta[,1]),6]
ord2 <- paste0(v1,':', v2)

expr1 = expr[, ord1]
design <- readRDS(here('covid','data','tex','design.rds'))
cellanno <- data.frame(cell = colnames(expr1), sample = sub(':.*', '', colnames(expr1)), stringsAsFactors = FALSE)
pt <- seq(1, length(ord1))
names(pt) <- ord1
saveRDS(expr1, here('covid','data','GSE155673_pbmc','useDiffusionMap','tex','log2norm.rds'))
saveRDS(pt, here('covid','data','GSE155673_pbmc','useDiffusionMap','tex','pseudotime.rds'))
saveRDS(cellanno, here('covid','data','GSE155673_pbmc','useDiffusionMap','tex','cellanno.rds'))
saveRDS(design, here('covid','data','GSE155673_pbmc','useDiffusionMap','tex','design.rds'))

expr2 = expr[, ord2]
design <- readRDS(here('covid','data','temra','design.rds'))
cellanno <- data.frame(cell = colnames(expr2), sample = sub(':.*', '', colnames(expr2)), stringsAsFactors = FALSE)
pt <- seq(1, length(ord2))
names(pt) <- ord2
saveRDS(expr2, here('covid','data','GSE155673_pbmc','useDiffusionMap','temra','log2norm.rds'))
saveRDS(pt, here('covid','data','GSE155673_pbmc','useDiffusionMap','temra','pseudotime.rds'))
saveRDS(cellanno, here('covid','data','GSE155673_pbmc','useDiffusionMap','temra','cellanno.rds'))
saveRDS(design, here('covid','data','GSE155673_pbmc','useDiffusionMap','temra','design.rds'))

