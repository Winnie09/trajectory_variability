library(data.table)
a = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/raw/GSE158769_meta_data.txt', header = T, data.table = F)

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/')
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex_addBatch/'

design = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')
a = a[a$donor %in% rownames(design), ]

## make a dataframe to store the sample batch, sex info
d = data.frame(batch = a[,5], donor = a[,11], stringsAsFactors = FALSE)
d = d[!duplicated(d[,2]), ]
sex = rep('male', nrow(d))
sex[design[d[,2], 2]==1] = 'female'
d = cbind(d, sex)
d[,3] = as.character(d[,3])

## only retain batches that have both sexes in the batch (should be at least 1 sample to contribute to the estimation of batch effects, should be at least 1 male and 1 female to contribute to the estimation of sex difference)
tb = table(d[,c(1,3)])
sum(tb[,1] > 0 & tb[,2] > 0)
tb.sel = tb[tb[,1] > 0 & tb[,2] > 0,]
batch.sel = rownames(tb.sel)
rownames(tb.sel) = paste0('batch', rownames(tb.sel))

write.csv(tb.sel, paste0(rdir, 'sample_batch_sex.csv'))

## update the dataframe of sample batch, sex info
d.select = d[d[,1] %in% batch.sel, ]
row.names(d.select) = d.select[,2]

## check how many cells are retained
cellanno = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
cellanno.select = cellanno[cellanno[,2] %in% d.select[,2], ]


##### =======================
##### prepare new data
##### =======================
batch = as.character(d.select[,1])
m = model.matrix(~batch)
rownames(m) = rownames(d.select)

design.sel = design[d.select[,2], ]
m = cbind(m, sex = design.sel[,2])
design2 = m

## paste batch to sample names in design
rownames(design2) = paste0(d.select[,2], '_batch', d.select[,1])
saveRDS(design2, paste0(rdir, 'design.rds'))

## paste batch to sample names in cellanno
v = paste0(cellanno.select[,2], '_batch', d.select[cellanno.select[,2],1])
cellanno2 = cellanno.select
cellanno2[,2] = v
saveRDS(cellanno2, paste0(rdir, 'cellanno.rds'))

## check
identical(sort(unique(v)), sort(rownames(design2)))
           
## expr
expr = readRDS('expr.rds')
expr2 = expr[, cellanno2[,1]]
saveRDS(expr2, paste0(rdir, 'expr.rds'))

## pt
ptpc2 = readRDS('ptpc2.rds')
ptpc4 = readRDS('ptpc4.rds')
saveRDS(ptpc2[cellanno2[,1]], paste0(rdir, 'ptpc2.rds'))
saveRDS(ptpc4[cellanno2[,1]], paste0(rdir, 'ptpc4.rds'))


