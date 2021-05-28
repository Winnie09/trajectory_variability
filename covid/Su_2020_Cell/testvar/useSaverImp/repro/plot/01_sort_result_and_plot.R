compareid <- (1:50)*100
scompareid <- (1:2)*100

af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/res',full.names = T)
d1 <- lapply(af,readRDS)
expid <- expand.grid(1:length(d1),1:length(d1))
expid <- expid[expid[,1] < expid[,2],]
ao1 <- sapply(1:nrow(expid),function(i) {
  s1 <- data.frame(d1[[unlist(expid[i,1])]][[1]])
  s2 <- data.frame(d1[[unlist(expid[i,2])]][[1]])
  s1 <- s1[order(s1$pval.overall,-s1$z.overall),]
  s2 <- s2[order(s2$pval.overall,-s2$z.overall),]
  mean(sapply(scompareid,function(i) {
    length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
  }))
}) 
i <- which.max(ao1)
s1 <- data.frame(d1[[unlist(expid[i,1])]][[1]])
s2 <- data.frame(d1[[unlist(expid[i,2])]][[1]])
s1 <- s1[order(s1$pval.overall,-s1$z.overall),]
s2 <- s2[order(s2$pval.overall,-s2$z.overall),]
o1 <- sapply(compareid,function(i) {
  length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
})

af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/tradeseq',full.names = T)
d2 <- lapply(af,readRDS)
expid <- expand.grid(1:length(d2),1:length(d2))
expid <- expid[expid[,1] < expid[,2],]
ao2 <- sapply(1:3,function(metid) {
  sapply(1:nrow(expid),function(i) {
    s1 <- data.frame(d2[[unlist(expid[i,1])]][[metid]])
    s2 <- data.frame(d2[[unlist(expid[i,2])]][[metid]])
    mean(sapply(scompareid,function(i) {
      length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
    }))
  }) 
})

i <- which.min(apply(ao2,1,max))

o2 <- sapply(1:3,function(metid) {
  s1 <- data.frame(d2[[unlist(expid[i,1])]][[metid]])
  s2 <- data.frame(d2[[unlist(expid[i,2])]][[metid]])
  o <- sapply(compareid,function(i) {
    length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
  })
})


af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/limma',full.names = T)[1:30]
ao3 <- sapply(af,function(i) {
  res <- readRDS(i)
  s1 <- res[[1]]
  s2 <- res[[2]]
  mean(sapply(scompareid,function(i) {
    length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
  }))
}) 

i <- af[which.min(ao3)]
res <- readRDS(i)
s1 <- res[[1]]
s2 <- res[[2]]
o3 <- sapply(compareid,function(i) {
  length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
})

df <- cbind(o1,o2,o3)
colnames(df) <- c('Lamian',paste0('tradeseq_',names(d2[[1]])),'Limma')
df <- melt(df)
colnames(df) <- c('numgene','method','overlap')
df[,1] <- compareid
saveRDS(df,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/df.rds')


df = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/df.rds')
library(reshape2)
library(ggplot2)
library(RColorBrewer)
mycol = brewer.pal(8, 'Dark2')[1:length(unique(df[,2]))]
mycol = mycol[c(2,1,3:5)]
names(mycol) = unique(as.character(df[,2]))
pdf('/Users/wenpinhou/Dropbox/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/overlapline.pdf',width=4,height=2.1)
ggplot(df,aes(x=numgene,y=overlap,col=method)) + 
  geom_point(size = 0.2) + 
  geom_line(size = 0.1) + theme_classic() +
  xlab('number of top genes') + ylab('overlap proportion') +
  scale_color_manual(values = mycol) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pd <- tapply(df[,3],df[,2],mean)
pd <- data.frame(method=names(pd),overlap=pd)
pdf('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/plot/overlap.pdf',width=5,height=5)
ggplot(pd,aes(method,overlap)) + geom_bar(stat='identity') + theme_classic()
dev.off()

