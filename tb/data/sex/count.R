library(data.table)
d <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/raw/GSE158769_exprs_raw.tsv',data.table=F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])

m <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/raw/GSE158769_meta_data.txt',data.table=F)
identical(sort(m[,1]),sort(colnames(d)))
colnames(d) <- paste0(m[match(colnames(d),m[,1]),'donor'],':',colnames(d))

e <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/expr.rds')

d <- d[rownames(e),colnames(e)]

saveRDS(d,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/count.rds')


