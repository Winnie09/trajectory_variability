library(ggplot2)
library(RColorBrewer)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/xci/plot/'

l <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/xci/list.csv',header=F)

## read in results
res <- do.call(rbind,sapply(c('erythroid','monocyte','lymph'),function(tis) {
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/',tis,'/gender/gender_res.rds'))[[1]]
  g <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/EM_pm/',tis,'/gender/differential_genes.csv'),as.is=T,row.names = 1)
  g <- sub(':.*','',rownames(g))
  ng <- setdiff(sub(':.*','',rownames(d)),g)
  t1 <- l[match(intersect(g,l[,1]),l[,1]),2]
  t2 <- l[match(intersect(ng,l[,1]),l[,1]),2]
  data.frame(tissue=tis,sig=rep(c('significant','not significant'),each=4),type=rep(c('E','mostly E','S','mostly S'),2),num=c(sum(t1=='E'),sum(t1=='Mostly E'),sum(t1=='S'),sum(t1=='Mostly S'),sum(t2=='E'),sum(t2=='Mostly E'),sum(t2=='S'),sum(t2=='Mostly S')))  
},simplify = F))

## plot results
pdf(paste0(pdir, 'sig_nosig_gene_proportion_vs_xci_barplot.pdf'), width = 3.5, height = 1.8)
ggplot(res,aes(x=type,y=num,fill=sig)) + geom_bar(stat='identity',position='stack') + facet_wrap(~tissue) + coord_flip() + ylab('mumber of genes') + xlab('X-chromosome inactivation') + scale_fill_brewer(palette = 'Dark2')  + scale_y_continuous(breaks=c(0, 100, 200))
dev.off()


prop <- res[res$sig=='significant',]
prop[,4] <- res[res$sig=='significant',4]/(res[res$sig=='significant',4]+res[res$sig=='not significant',4])
pdf(paste0(pdir, 'xci_gene_proportion_barplot.pdf'), width = 2.8, height = 1.8)
ggplot(prop,aes(x=type,y=num)) + geom_bar(stat='identity',position='stack') + facet_wrap(~tissue) + coord_flip() + ylab('percentage of significant genes') + xlab('X-chromosome inactivation') + scale_y_continuous(breaks=c(0, 0.3, 0.5))
dev.off()


