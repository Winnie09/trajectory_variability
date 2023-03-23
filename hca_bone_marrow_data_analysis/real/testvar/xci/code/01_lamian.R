l <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/xci/list.csv',header=F)

res <- do.call(rbind,sapply(c('erythroid','monocyte','lymph'),function(tis) {
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/',tis,'/gender/gender_res.rds'))[[1]]
  g <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/EM_pm/',tis,'/gender/differential_genes.csv'),as.is=T,row.names = 1)
  g <- sub(':.*','',rownames(g))
  ng <- setdiff(sub(':.*','',rownames(d)),g)
  t1 <- l[match(intersect(g,l[,1]),l[,1]),2]
  t2 <- l[match(intersect(ng,l[,1]),l[,1]),2]
  data.frame(tissue=tis,sig=rep(c('sig','nosig'),each=4),type=rep(c('E','mostly E','S','mostly S'),2),num=c(sum(t1=='E'),sum(t1=='Mostly E'),sum(t1=='S'),sum(t1=='Mostly S'),sum(t2=='E'),sum(t2=='Mostly E'),sum(t2=='S'),sum(t2=='Mostly S')))  
},simplify = F))

library(ggplot2)
ggplot(res,aes(x=type,y=num,fill=sig)) + geom_bar(stat='identity',position='stack') + facet_wrap(~tissue) + coord_flip() + theme_classic()

prop <- res[res$sig=='sig',]
prop[,4] <- res[res$sig=='sig',4]/(res[res$sig=='sig',4]+res[res$sig=='nosig',4])
ggplot(prop,aes(x=type,y=num)) + geom_bar(stat='identity',position='stack') + facet_wrap(~tissue) + coord_flip() + theme_classic()
