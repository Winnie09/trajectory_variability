d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/proc/pca.rds')
clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/proc/cluster.rds')
ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/ct/sc.rds')
ct <- ct[rownames(d)]
library(ggplot2)
library(gridExtra)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hcahar/data/plot/ct.pdf',width=12,height=4)
p1 <- ggplot(data.frame(pc1=d[,1],pc2=d[,2],ct=ct),aes(x=pc1,y=pc2,col=ct)) + geom_point(size=0.2) + theme_classic()
p2 <- ggplot(data.frame(pc1=d[,1],pc3=d[,3],ct=ct),aes(x=pc1,y=pc3,col=ct)) + geom_point(size=0.2) + theme_classic()
p3 <- ggplot(data.frame(pc2=d[,2],pc3=d[,3],ct=ct),aes(x=pc2,y=pc3,col=ct)) + geom_point(size=0.2) + theme_classic()
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

