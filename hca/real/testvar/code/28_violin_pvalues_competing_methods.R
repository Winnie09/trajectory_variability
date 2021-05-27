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

for (celltype in setdiff(list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm'),'perf')) {
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
  
  sig <- list(Lamian_excludeLimma = setdiff(lamiansig,limmasig),tradeSeq_excludeLamian = setdiff(unlist(tradeseqsig),lamiansig))
  
  allg <- sub(':.*','',rownames(tradeseq[[1]]))
  
  permud <- reald <- NULL
  for (met in names(sig)) {
    gn <- sig[[met]]
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
  
  for (met in names(genes)) {
    gn <- genes[[met]]
    allg <- sub(':.*','',gn)
    v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
    v1 <- mean(v1[seq(1,length(v1)) %% 10 == 0])
    v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
    v2 <- mean(v2[seq(1,length(v2)) %% 10 == 0])
    # saveRDS(v1, paste0(rdir, path, '/gender/gender_chrX_overlap.rds'))
    # saveRDS(v2, paste0(rdir, path, '/gender/gender_chrY_overlap.rds'))
    
    ##### permute reported gene order
    v1_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
      set.seed(myseed+100)
      allg.pm = sample(allg)
      tmp = cumsum(allg.pm %in% u1)/seq(1,length(allg.pm))
      mean(tmp[seq(1,length(tmp)) %% 10 == 0])
    },mc.cores=detectCores()))
    
    
    v2_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
      set.seed(myseed+100)
      allg.pm = sample(allg)
      tmp = cumsum(allg.pm %in% u2)/seq(1,length(allg.pm))
      mean(tmp[seq(1,length(tmp)) %% 10 == 0])
    },mc.cores=detectCores()))
    
    permud <- rbind(permud,data.frame(per = c(v1_pm,v2_pm), type=rep(c('chrX','chrY'),each=1e4),method=met, stringsAsFactors = FALSE))
    reald <- rbind(reald,data.frame(per=c(v1,v2),type=c('chrX','chrY'),pvalue=c(mean(v1_pm >= v1),mean(v2_pm >= v2)),method=met,stringsAsFactors = F))
  }
  permud$method <- factor(permud$method,levels=c('Lamian_excludeLimma', 'tradeSeq_excludeLamian', 'Lamian','limma', 'tradeSeq_diffEndTest', 'tradeSeq_patternTest', 'tradeSeq_earlyDETest'))
  reald$method <- factor(reald$method,levels=c('Lamian_excludeLimma', 'tradeSeq_excludeLamian', 'Lamian','limma', 'tradeSeq_diffEndTest', 'tradeSeq_patternTest', 'tradeSeq_earlyDETest'))
  
  
  library(RColorBrewer)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin.pdf'),width=5,height=5)
  print(ggplot() + geom_violin(data=permud,aes(x=method,y=per,col=type)) + geom_point(data=reald,aes(x=method,y=per,col=type),size=3) + geom_text(data=reald,aes(x=method,y=max(reald$per)*1.3,label=pvalue)) + theme_classic() + facet_wrap(~type) + coord_flip(ylim=c(0,max(reald$per)*1.5)) + xlab('') + ylab('Proportion') + scale_color_manual(values=c('chrX'=brewer.pal(3,'Pastel1')[1],'chrY'=brewer.pal(3,'Pastel1')[2])) + theme(legend.position = 'bottom',strip.background = element_blank(),strip.text = element_text(size=20),legend.title = element_blank()) + scale_y_continuous(breaks=c(0,round(max(reald$per)*0.6,2),round(max(reald$per)*1.2,2))))
  dev.off()
  
}

