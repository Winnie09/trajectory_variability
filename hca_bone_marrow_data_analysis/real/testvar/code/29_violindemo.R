rm(list=ls())
library(here)
library(ggplot2)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

celltype <- 'monocyte'
lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/',celltype,'/gender/gender_res.rds'))[[1]]
lamiangene <- rownames(lamian)[order(lamian[,2])]
if (celltype=='monocyte') {
  lamiansig <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]  
} else {
  lamiansig <- rownames(lamian)[lamian[,'pval.overall'] < 0.05]    
}
limma <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/limma/',celltype,'/gender_res.rds'))
limmagene <- rownames(limma)[order(limma$P.Value)]
limmasig <- rownames(limma)[limma$adj.P.Val < 0.05]
tradeseq <- readRDS(paste0('hca/real/testvar/result/tradeSeq/',celltype,'/gender/testvar_res.rds'))
tradeseqgene <- lapply(tradeseq,function(i) rownames(i)[order(i$P.Value,-i$waldStat)])
tradeseqsig <- lapply(tradeseq,function(i) rownames(i)[i$adj.P.Val < 0.05])

genes <- list(Lamian=lamiangene)
for (i in names(tradeseq))
  genes[[paste0('tradeSeq_',i)]]=tradeseqgene[[i]]
genes[['limma']] <- limmagene

pd <- NULL
for (met in names(genes)) {
  gn <- genes[[met]]
  allg <- sub(':.*','',gn)
  pd <- rbind(pd,data.frame(ng=cumsum(allg %in% u1)/seq(1,length(allg))*100,tg=seq(1,length(allg)),met=met,stringsAsFactors = F))
}

set.seed(1234)
gn <- sample(genes[[1]])
allg <- sub(':.*','',gn)
pd <- rbind(pd,data.frame(ng=cumsum(allg %in% u1)/seq(1,length(allg))*100,tg=seq(1,length(allg)),met='random',stringsAsFactors = F))

cv <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/palette/cv.rds')
cv <- c(cv,'random'='grey')

library(ggplot2)
library(RColorBrewer)
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/violindemo_allnum.pdf'),width=4,height=3)
ggplot(pd[pd$tg <= 1000,],aes(x=tg,y=ng,col=met,group=met)) + geom_line() + theme_classic() + scale_color_manual(values=cv) + theme(legend.title = element_blank()) + xlab('number of top genes') + ylab('proportion of X chromosome genes (%)')
dev.off()


pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/violindemo.pdf'),width=4,height=3)
ggplot(pd[pd$tg <= 1000 & pd$tg %% 100==0,],aes(x=tg,y=ng,col=met,group=met)) + geom_line() + theme_classic() + scale_color_manual(values=cv) + theme(legend.title = element_blank()) + xlab('number of top genes') + ylab('overlap proportion (%)')
dev.off()



