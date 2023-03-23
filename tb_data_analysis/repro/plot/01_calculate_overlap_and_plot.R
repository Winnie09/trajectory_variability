library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(here)
library(parallel)
setwd(here())
compareid <- (1:50)*100
scompareid <- (1:2)*100
rdir <- 'tb/repro/summary/'

res <- readRDS(paste0(rdir, 'Lamian.rds'))
Lamian.pm.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

res <- readRDS(paste0(rdir, 'Lamian.chisq.rds'))
Lamian.chisq.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

res <- readRDS(paste0(rdir, 'limma.rds'))
limma.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

### monocle2 trajtest correction
af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/monocle2_trajtest_corr',full.names = T)
ao5 <- sapply(1:5,function(i) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/monocle2_trajtest_corr/',i, '_1.rds'))
  s1 = s1[order(s1[,3], -s1[, 1]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/monocle2_trajtest_corr/',i, '_2.rds'))
  s2 = s2[order(s2[,3], -s2[, 1]),]
  mean(sapply(scompareid,function(j) {
    length(intersect(rownames(s1)[1:j],rownames(s2)[1:j]))/j
  }))
})
i <- which.min(ao5)
s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/monocle2_trajtest_corr/',i, '_1.rds'))
s1 = s1[order(s1[,3], -s1[, 1]),]
s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/monocle2_trajtest_corr/',i, '_2.rds'))
s2 = s2[order(s2[,3], -s2[, 1]),]
saveRDS(list(s1, s2), paste0(rdir, 'monocle2_trajTest.corr.rds'))
res <- readRDS(paste0(rdir, 'monocle2_trajTest.corr.rds'))
monocle2_trajTest.corr.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

saveRDS(list(Lamian.pm=Lamian.pm.o,Lamian.chisq = Lamian.chisq.o, limma=limma.o, monocle2_trajTest.corr=monocle2_trajTest.corr.o), paste0(rdir, 'overlap.rds'))

## read in plot data
o = readRDS(paste0(rdir, 'overlap.rds'))
df = do.call(cbind, o)
colnames(df)[4] <- 'monocle2TrajtestCorr'

coldf = read.csv('resource/color_code.csv', as.is = T, row.names = 1)
colv = coldf[,2]
names(colv) = coldf[,3]

df <- reshape2::melt(df)
colnames(df) <- c('numgene','method','overlap')
df[,1] <- compareid
saveRDS(df,paste0(rdir,'/df.rds'))

##########################################################################3
df = readRDS(paste0(rdir,'/df.rds'))
source('resource/ggplot_theme.R')
theme_set(.new_theme)

pdf('tb/repro/plot/overlapline.pdf',width=4,height=2.2)
ggplot(df,aes(x=numgene,y=overlap,col=method, shape = method)) + 
  geom_point(size = 0.2) + 
  geom_line(size = 0.1) + theme_classic() +
  xlab('number of top genes') + ylab('overlap proportion') +
  scale_color_manual(values = colv) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pd <- tapply(df[,3],df[,2],mean)
pd <- data.frame(method=names(pd),overlap=pd)
pdf('tb/repro/plot/avecov.pdf',width=1.5,height=2)
ggplot(pd,aes(method,overlap, fill = method), alpha = 0.1) + geom_bar(stat='identity') + 
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), legend.position = 'none') +
  scale_color_manual(values = colv) +
  ylab('average overlap') +
  xlab('') 
dev.off()


