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
tradeseq <- readRDS(paste0('hca/real/testvar/result/tradeSeq/monocyte/gender/testvar_res.rds'))

for (sn in names(tradeseq)) {
  res <- tradeseq[[sn]]
  gn <- setdiff(rownames(res)[res$adj.P.Val < 0.05],lamian)
  gn <- sub(':.*','',gn)
  v1 <- mean(gn %in% u1)
  v2 <- mean(gn %in% u2)
  
  allg <- sub(':.*','',rownames(res))
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
  
  pd <- data.frame(chrX_pm = v1_pm, chrY_pm = v2_pm, stringsAsFactors = FALSE)
  p1 <- ggplot(data = pd, aes(x=chrX_pm)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") +
    geom_vline(xintercept = v1, color = 'red') +
    theme_classic()+
    xlab('Overlap proportion mean')+ ylab('Density')+
    ggtitle(paste0('chromosome X, p=', round(mean(v1_pm > v1),3)))
  
  p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
    geom_density(alpha=.2, fill="lightblue") +
    geom_vline(xintercept = v2, color = 'darkblue') +
    theme_classic()+
    xlab('Overlap proportion mean')+ ylab('Density') + 
    ggtitle(paste0('chromosome Y, p=', round(mean(v2_pm >v2),3)))
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/tradeSeq_excEM/',sn,'.pdf'), width = 5, height = 2.1)
  gridExtra::grid.arrange(p1,p2,nrow=1)
  dev.off()  
}


