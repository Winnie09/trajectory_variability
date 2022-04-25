library(Seurat)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcarev/data/count.rds')
clu <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcarev/data/harmony/cluster.rds'))
clu <- clu[colnames(d)]
dm <- sapply(1:max(clu),function(i) rowSums(d[,clu==i]))
rc <- colSums(dm)
rc <- rc/median(rc)
dm <- log2(t(t(dm)/rc) + 1)
row.names(dm) <- sub(':.*','',row.names(dm))
dm <- dm[!duplicated(row.names(dm)),]

data <- read.table("/home-4/zji4@jhu.edu/scratch/scdata/HCA/geobulk/GSE74246/data/GSE74246_RNAseq_All_Counts.txt",header=T)
row.names(data) <- data[,1]
data <- as.matrix(data[,-1])
data <- data[,grep("^X",colnames(data))]
gl <- readRDS("/home-4/zji4@jhu.edu/scratch/resource/gl/res/hg19.rds")
gl <- gl/1000
load("/home-4/zji4@jhu.edu/scratch/resource/gn/res/hg19.rda")
tgngl <- tapply(gl[geneid],genename,max)
gngl <- as.vector(tgngl)
names(gngl) <- names(tgngl)
data <- data[row.names(data) %in% names(gngl),]
data <- data/gngl[row.names(data)]

libsize <- colSums(data)
libsize <- libsize/1e6
data <- sweep(data,2,libsize,"/")
data <- log2(data + 1)
data <- data[rowSums(data >= 1) >=1,]

ct <- sub('.*\\.','',colnames(data))
bd <- sapply(unique(ct),function(sct) {
  rowMeans(data[,ct==sct])
})


source('/home-4/zji4@jhu.edu/scratch/raisin/software/proc/ct.R')
ct <- ctfunc()
saveRDS(ct,file='/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcarev/data/harmony/ct.rds')


