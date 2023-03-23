library(parallel)
library(limma)

sg <- row.names(readRDS('/home-4/zji4@jhu.edu/scratch/raisin/HCA/common/proc/matrix/count.rds'))
sg <- sub(':.*','',sg)
data <- read.table("/home-4/zji4@jhu.edu/scratch/scdata/HCA/geobulk/GSE74246/data/GSE74246_RNAseq_All_Counts.txt",header=T)
row.names(data) <- data[,1]
data <- as.matrix(data[,-1])
data <- data[,grep("^X",colnames(data))]
ct <- sub('.*\\.','',colnames(data))
alldata <- data

alldata <- alldata[,sub('.*\\.','',colnames(alldata)) %in% c('HSC','Mono','CD4Tcell','CD8Tcell','NKcell','Bcell','Ery')]
alldata <- alldata[rowSums(alldata) > 0,]
alldata <- t(log2(t(alldata)/colSums(alldata) * 1e6 + 1))

ct <- sub('.*\\.','',colnames(alldata))
ct[ct %in% c('CD4Tcell','CD8Tcell','NKcell')] <- 'TNK'
gl <- sapply(unique(ct),function(sct) {
  des <- cbind(1,ct==sct)
  res <- topTable(eBayes(lmFit(alldata,des)),coef=2,n=nrow(alldata))
  res <- res[res[,1] > 0,]  
  res <- res[order(-res[,1]),]
  rownames(res)[1:100]
},simplify = F)
gl[['Lym']] <- union(gl[['TNK']],gl[['Bcell']])
gl <- gl[!names(gl) %in% c('TNK','Bcell')]

saveRDS(gl, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/bulkdiff/result/diffgene.rds')






