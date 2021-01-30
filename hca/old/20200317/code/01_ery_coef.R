setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/data/HCA/')
library(Seurat)
set.seed(12345)
u = readRDS('./proc/integrate/umap/umap.rds')
clu = readRDS('./proc/cluster/resolution0.1.rds')
n <- names(clu)
clu = as.numeric(as.character(clu))
names(clu) = n

library(TSCAN)
set.seed(12345)
mcclu <- as.numeric(as.factor(clu[clu %in% c(7,14,8,3,0,4)]))
names(mcclu) <- row.names(u[clu %in% c(7,14,8,3,0,4),])
mc <- exprmclust(t(u[clu %in% c(7,14,8,3,0,4),]),cluster=mcclu,reduce=F)

library(igraph)
set.seed(12345)
MSTdf = data.frame(from=c(4,6,5,2,1), to=c(6,5,2,1,3))
mc$MSTtree = graph_from_data_frame(MSTdf, directed = FALSE, vertices = NULL)
order <- TSCANorder(mc,orderonly = T,MSTorder=c(4,6,5,2,1,3))
order = data.frame(Pseudotime=1:length(order), Cell=order,stringsAsFactors = F)
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/ery/order.rds')

mat = readRDS('./proc/matrix/saver.rds')
mat = mat[,order$Cell]
meta = read.table('./meta/meta.tsv',sep='\t',header = T)
meta = meta[grepl('matrix',meta[,'analysis_file.file_core.format']),]
v = meta[,'donor_organism.sex']
v = sapply(rownames(u), function(i) sub(':.*','',i))
v = gsub('BM','Bone Marrow donor ',v)
v = v[match(colnames(mat), names(v))]

res <- sapply(unique(v), function(s){
  print(s)  
  trainData = mat[, v[match(colnames(mat),names(v))]==s]
  trainX = order[match(colnames(trainData), order$Cell),'Pseudotime']
  
  num.base = 10
  knots = seq(min(order$Pseudotime),max(order$Pseudotime),length.out=num.base+2)[2:(num.base+1)]
  library(splines)
  base = cbind(1,bs(trainX,knots = knots))
  colidx = NULL
  for (ii in seq(2,ncol(base))){
    if (length(unique(base[,ii])) == 1) colidx = c(colidx, ii)
  }
  if (length(colidx)) base = base[,-colidx]
  para = chol2inv(chol(crossprod(base))) %*% t(base) %*% t(trainData) ##
  saveRDS(t(para), paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/result/HCA/ery/spline_coef_', gsub(':','_',sub('_.*','',colnames(trainData)[1])), '.rds'))
  return(0)
})

# ## check para  
# lm(trainData[1,]~base)
# ## check marker genes
# gid  = grep('^HBA2:', rownames(mat))
# plot(trainData[gid,]~trainX)
# y = t(base %*% para)[gid,]
# lines(y~trainX, col='red')
