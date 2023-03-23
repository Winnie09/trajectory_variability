rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
pdir <- 'tb/plot/sex/'

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home/whou10/scratch16/whou10/resource/chrX_genename.rds')
u2 = readRDS('/home/whou10/scratch16/whou10/resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')
l <- read.csv('hca/real/xci/list.csv',header=F)

## Lamian.pm
# lamiangene <- readRDS('/home/whou10/scratch16/whou10/trajectory_variability/tb/res/sex/pc2/lamian_pm_allcores.rds')[[1]]
# lamiangene <- rownames(lamiangene)
# lamiansig <- rownames(read.csv(paste0(pdir, 'differential_genes.csv'),as.is=T,row.names = 1))
# 
# ## Lamian.chisq
# lamc <- readRDS('/home/whou10/scratch16/whou10/trajectory_variability/tb/res/sex/pc2/lamian_chisq.rds')[[1]]
# lamcgenes <- rownames(lamc)[order(lamc[,2])] ## all genes order
# lamcsig <- rownames(lamc)[lamc[,'fdr.chisq.overall'] < 0.05]   ## significant genes

# genes <- list(Lamian.pm=lamiangene)
# genes[['Lamian.chisq']] <- lamcgenes
# 
# sig <- list(Lamian.pm=lamiansig)
# sig[['Lamian.chisq']] <- lamcsig


sig <- readRDS(paste0('/home/whou10/scratch16/whou10/trajectory_variability/tb/perf/violin_plotdata_genesSig.rds'))
names(sig)[1] = 'Lamian.pm'

genes <- readRDS(paste0('/home/whou10/scratch16/whou10/trajectory_variability/tb/perf/violin_plotdata_genes.rds'))
names(genes)[which(names(genes) == 'Lamian')] <- 'Lamian.pm'
names(genes)[which(names(genes) == 'monocle2trajtest')] <- 'monocle2Trajtest'
names(genes)[which(names(genes) == 'monocle2trajtest.corr')] <- 'monocle2TrajtestCorr'
genes <- genes[1:5]


res <- sapply(names(sig),function(n) {
  g <- sub(':.*','',sig[[n]])
  ng <- setdiff(sub(':.*','',genes[[n]]),g)
  t1 <- l[match(intersect(g,l[,1]),l[,1]),2]
  t2 <- l[match(intersect(ng,l[,1]),l[,1]),2]
  mean(t1 %in% c('E','Mostly E'))
})
res = rev(sort(res))
pd = data.frame(method = names(res), value = res, stringsAsFactors = F)
pd[,1] = factor(pd[,1], levels = names(res))
saveRDS(pd, paste0(pdir, 'xci_compare_pd.rds'))

source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
# source('/Users/wenpinhou/Dropbox/trajectory_variability/resource/ggplot_theme.R')

library(ggplot2)
pdf(paste0(pdir, '/xci_compare.pdf'), height = 1.5, width = 3.5)
ggplot(pd,aes(x=method,y=value)) + geom_bar(stat='identity')  + coord_flip() + ylab('Percentage of differential genes escaping X-inactivation') + xlab('')
dev.off()

##############################################
## read in results
d <- genes[['Lamian.pm']]
g <- sig[['Lamian.pm']]
ng <- setdiff(d,g)
t1 <- l[match(intersect(g,l[,1]),l[,1]),2]
t2 <- l[match(intersect(ng,l[,1]),l[,1]),2]
res <- data.frame(sig=rep(c('significant','not significant'),each=4),type=rep(c('E','mostly E','S','mostly S'),2),num=c(sum(t1=='E'),sum(t1=='Mostly E'),sum(t1=='S'),sum(t1=='Mostly S'),sum(t2=='E'),sum(t2=='Mostly E'),sum(t2=='S'),sum(t2=='Mostly S')))  
saveRDS(res, paste0(pdir, 'xci_sig_nosig_gene_proportion_vs_xci_barplot_pd.rds'))

## plot results
pdf(paste0(pdir, 'xci_sig_nosig_gene_proportion_vs_xci_barplot.pdf'), width = 3.5, height = 1.8)
ggplot(res,aes(x=type,y=num,fill=sig)) + geom_bar(stat='identity',position='stack')  + coord_flip() + ylab('number of genes') + xlab('X-chromosome inactivation') + scale_fill_brewer(palette = 'Dark2')  + scale_y_continuous(breaks=c(0, 100, 200))
dev.off()

prop <- res[res$sig=='significant',]
prop[,3] <- res[res$sig=='significant',3]/(res[res$sig=='significant',3]+res[res$sig=='not significant',3])
saveRDS(prop, paste0(pdir, 'xci_gene_proportion_barplot_pd.rds'))

pdf(paste0(pdir, 'xci_gene_proportion_barplot.pdf'), width = 3, height = 1.5)
ggplot(prop,aes(x=type,y=num)) + geom_bar(stat='identity',position='stack') + coord_flip() + ylab('Percentage of significant genes') + xlab('X-chromosome inactivation') + scale_y_continuous(breaks=c(0, 0.3, 0.5, 0.7))
dev.off()

