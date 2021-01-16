library(here)
library(ggplot2)
setwd(here())
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')

for (path in c('erythroid', 'lymph', 'monocyte')){
  print(path)
  res <- readRDS(paste0(rdir, path, '/gender_fdr_res.rds'))
  allg <- sub(':.*', '', rownames(res))
  res <- res[res[,1] < 0.05, ]
  diffgene <- sub(':.*', '', rownames(res))
  res.lm <- readRDS(paste0(ddir, path,'/meandiff_gender_res.rds'))
  res.lm <- res.lm[res.lm[,5] < 0.05, ]
  diffgene.lm <- sub(':.*', '', rownames(res.lm))
  
  newgene <- setdiff(diffgene, diffgene.lm)
  v1 <- mean(newgene %in% u1 )
  v2 <- mean(newgene %in% u2 )
  v1_pm <- sapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(allg, length(intersect(newgene,u1)))
    mean(diffgene %in% w1)
  })
  
  v2_pm <- sapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(allg, length(intersect(newgene,u2)))
    mean(diffgene %in% w1)
  })
  
  x.pval <- mean(v1_pm > v1)
  y.pval <- mean(v2_pm > v2)
  
  pd <- data.frame(chrX_pm = v1_pm, chrY_pm = v2_pm, stringsAsFactors = FALSE)
  p1 <- ggplot(data = pd, aes(x=chrX_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white")+
   geom_density(alpha=.2, fill="#FF6666") +
   geom_vline(xintercept = v1, color = 'red') +
    theme_classic()+
    xlab('newgene overlap proportion mean')+
    ggtitle(paste0('chrX, p.value = ', x.pval))
  
  p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
   geom_density(alpha=.2, fill="lightblue") +
   geom_vline(xintercept = v1, color = 'darkblue') +
    theme_classic()+
    xlab('newgene overlap proportion mean')+
    ggtitle(paste0('chrY, p.value = ', y.pval))
  
  pdf(paste0(pdir, path, '/gender_hist_newgene_overlap_proportion_mean.pdf'), width = 7, height = 3)
  gridExtra::grid.arrange(p1,p2,nrow=1)
  dev.off()
}

rm(list=ls())
