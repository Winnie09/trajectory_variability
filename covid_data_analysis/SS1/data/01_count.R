library(Matrix)
m <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/data/EGAS00001004571_pbmc/meta/trans.rds')
d <- readRDS('/home-4/zji4@jhu.edu/scratch/covid/data/EGAS00001004571_pbmc/proc/10x.rds')
p1 <- sub(':.*','',colnames(d))
p2 <- sub('.*:','',colnames(d))
print(mean(p1 %in% m[,1]))
colnames(d) <- paste0(m[match(p1,m[,1]),2],':',p2)


id <- 1:ncol(d)
gn <- colSums(d > 0)  
id <- intersect(id,which(gn >= 500))

rc <- colSums(d)
mito <- colSums(d[grep('^MT-',row.names(d),ignore.case=T),])/rc
id <- intersect(id,which(rc >= 2000))
id <- intersect(id,which(mito <= 0.2))

d <- d[,id]
p <- sub(':.*','',colnames(d))
tab <- table(p)
d <- d[,p %in% names(tab)[tab >= 500]]
d <- d[rowSums(d) > 0,]
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/diffpt/covid/ss/count.rds')

