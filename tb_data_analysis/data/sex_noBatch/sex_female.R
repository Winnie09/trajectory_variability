library(data.table)
expr <- fread(paste0('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/tb/data/raw/GSE158769_exprs_raw.tsv'),data.table=F)
m <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/raw/GSE158769_meta_data.txt',data.table=F)
rownames(expr) <- expr[,1]
expr <- as.matrix(expr[,-1])

int <- intersect(colnames(expr),m[,1])
expr <- expr[,int]
m <- m[match(int,m[,1]),]
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gx <- sub('\".*','',sub('.*gene_name \"','',gtf[gtf[,1]=='chrX',9]))
gy <- sub('\".*','',sub('.*gene_name \"','',gtf[gtf[,1]=='chrY',9]))
gy <- setdiff(gy,gx)

## chromosome Y gene expression: cc
cc <- sapply(unique(m$donor),function(i) {
  sum(expr[intersect(gy,rownames(expr)),m[m$donor==i,1]])/sum(expr[,m[m$donor==i,1]])
})

## if chromasome Y gene expression smaller than a cutoff, then it is a female sample
female <- cc < 0.0004
saveRDS(female,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex_female.rds')

