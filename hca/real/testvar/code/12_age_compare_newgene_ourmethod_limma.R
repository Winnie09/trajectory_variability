library(here)
library(ggplot2)
setwd(here())
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'
u1 <- readRDS('agediff/result/age_diffgene.rds')


for (path in c('erythroid', 'lymph', 'monocyte')){
  print(path)
  res <- readRDS(paste0(rdir, path, '/age_fdr_res.rds'))
  res.lm <- readRDS(paste0(ddir, path,'/meandiff_age_res.rds'))
  allg <- sub(':.*', '', rownames(res))
  res <- res[res[,1] < 0.05, ]
  print(str(res))
  diffgene <- sub(':.*', '', rownames(res))
  res.lm <- res.lm[res.lm[,5] < 0.05, ]
  print(str(res.lm))
  diffgene.lm <- sub(':.*', '', rownames(res.lm))
  
  newgene <- setdiff(diffgene, diffgene.lm)
  v1 <- mean(newgene %in% u1 )
  
  v1_pm <- sapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(allg, length(intersect(newgene,u1)))
    mean(diffgene %in% w1)
  })
  pval <- mean(v1_pm > v1)

  pd <- data.frame(chrX_pm = v1_pm, stringsAsFactors = FALSE)
  pdf(paste0(pdir, path, '/age_hist_newgene_overlap_proportion_mean.pdf'), width = 3.5, height = 3)
  print(ggplot(data = pd, aes(x=chrX_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white")+
   geom_density(alpha=.2, fill="#FF6666") +
   geom_vline(xintercept = v1, color = 'red') +
    theme_classic()+
    xlab('newgene overlap proportion mean')+
    ggtitle(paste0('age, p.value = ', pval)))
  dev.off()
}

rm(list=ls())
