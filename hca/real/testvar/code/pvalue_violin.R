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

lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/monocyte/gender/gender_res.rds'))[[1]]
lamian <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]
limma <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/limma/monocyte/gender_res.rds')
limma <- rownames(limma)[limma$adj.P.Val < 0.05]
tradeseq <- readRDS(paste0('hca/real/testvar/result/tradeSeq/monocyte/gender/testvar_res.rds'))

gs <- list(Lamian=lamian)
for (i in names(tradeseq))
  gs[[paste0('tradeSeq_',i)]]=rownames(tradeseq[[i]])[tradeseq[[i]]$adj.P.Val < 0.05]
gs[['Lamian_excludeLimma']] <- setdiff(lamian,limma)
gs[['tradeSeq_excludeLamian']] <- setdiff(unlist(gs[grep('tradeSeq_',names(gs))]),lamian)

allg <- sub(':.*','',rownames(tradeseq[[1]]))

permud <- reald <- NULL
for (met in names(gs)) {
  gn <- gs[[met]]
  gn <- sub(':.*','',gn)
  v1 <- mean(gn %in% u1)
  v2 <- mean(gn %in% u2)
  
  ##### permute reported gene order
  v1_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed+100)
    gnpm = sample(allg,length(gn))
    mean(gnpm %in% u1)
  },mc.cores=detectCores()))
  
  
  v2_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed+100)
    gnpm = sample(allg,length(gn))
    mean(gnpm %in% u2)
  },mc.cores=detectCores()))
  
  permud <- rbind(permud,data.frame(per = c(v1_pm,v2_pm), type=rep(c('chrX','chrY'),each=1e4),method=met, stringsAsFactors = FALSE))
  reald <- rbind(reald,data.frame(per=c(v1,v2),type=c('chrX','chrY'),pvalue=c(mean(v1_pm >= v1),mean(v2_pm >= v2)),method=met,stringsAsFactors = F))
}

library(RColorBrewer)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/pvalue_violin.pdf')
ggplot() + geom_violin(data=permud,aes(x=method,y=per,col=type)) + geom_point(data=reald,aes(x=method,y=per,col=type),size=3) + geom_text(data=reald,aes(x=method,y=max(reald$per)*1.25,label=pvalue)) + theme_classic() + facet_wrap(~type) + coord_flip(ylim=c(0,max(reald$per)*1.5)) + xlab('') + ylab('Proportion') + scale_color_manual(values=c('chrX'=brewer.pal(3,'Pastel1')[1],'chrY'=brewer.pal(3,'Pastel1')[2])) + theme(legend.position = 'bottom',strip.background = element_blank(),strip.text = element_text(size=20),legend.title = element_blank()) + scale_y_continuous(breaks=c(0,0.1,0.2))
dev.off()

