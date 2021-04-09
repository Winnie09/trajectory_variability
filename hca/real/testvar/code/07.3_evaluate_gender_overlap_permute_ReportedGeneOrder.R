rm(list=ls())
library(here)
library(ggplot2)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/EM_pm/'
pdir <- 'hca/real/testvar/plot/EM_pm/'

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

for (path in c('erythroid', 'lymph', 'monocyte')){
  ## read data
  Res <- readRDS(paste0(ddir, path, '/gender/gender_res.rds'))
  res <- Res$statistics
  ## identify DDGType 
  DDGType <- getDDGType(Res)
  # ------------------
  #  Evaluation Gender 
  # ------------------
  ## TP
  res <- res[order(res[,'fdr.overall'], -res[,'z.overall']), ]
  allg <- sub(':.*', '', rownames(res))
  diffgene = sub(':.*', '', rownames(res[res[,'fdr.overall'] < 0.05,]))
  str(diffgene)
  ### proportion with Sex gold standard
  #### using all DDGs
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
  
  names(v1_pm) <- paste0('top',seq(1,length(v1_pm)))
  names(v2_pm) <- paste0('top',seq(1,length(v2_pm)))
  
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
  
  pdf(paste0(pdir, path, '/gender/gender_hist_mean_of_overlap_proportion_PmReportedGeneOrder.pdf'), width = 5, height = 2.1)
  gridExtra::grid.arrange(p1,p2,nrow=1)
  dev.off()
  
  #### using only one type of DDG: trendSig, meanSig, other, bothSig
  if (length(diffgene) > 0){
    for (i in setdiff(unique(DDGType), 'nonDDG')){
      print(i)
      SigGene <- names(DDGType)[DDGType == i]
      if (i == 'meanSig'){
        allg <- sub(':.*','',intersect(rownames(res[order(res[,'fdr.meanDiff'], -res[,'z.meanDiff']), ]), SigGene))
      } else if (i == 'trendSig'){
        allg <- sub(':.*','',intersect(rownames(res[order(res[,'fdr.trendDiff'], -res[,'z.trendDiff']), ]), SigGene))
      } else if (i == 'bothSig' | i == 'other'){
        allg <- sub(':.*','',intersect(rownames(res[order(res[,'fdr.overall'], -res[,'z.overall']), ]), SigGene))
      }
      
      v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
      v1 <- mean(v1[seq(1,length(v1)) %% 10 == 0])
      v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
      v2 <- mean(v2[seq(1,length(v2)) %% 10 == 0])
      # saveRDS(v1, paste0(rdir, path, '/gender/gender_chrX_overlap_', i,'.rds'))
      # saveRDS(v2, paste0(rdir, path, '/gender/gender_chrY_overlap_', i, '.rds'))
      
      ##### permute gold standard genes
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
    
      names(v1_pm) <- paste0('top',seq(1,length(v1_pm)))
      names(v2_pm) <- paste0('top',seq(1,length(v2_pm)))
      
      pd <- data.frame(chrX_pm = v1_pm, chrY_pm = v2_pm, stringsAsFactors = FALSE)
      p1 <- ggplot(data = pd, aes(x=chrX_pm)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="#FF6666") +
        geom_vline(xintercept = v1, color = 'red') +
        theme_classic()+
        xlab('overlap proportion mean')+
        ggtitle(paste0('chromosome X, p=', round(mean(v1_pm > v1),3)))
      
      p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
        geom_density(alpha=.2, fill="lightblue") +
        geom_vline(xintercept = v2, color = 'darkblue') +
        theme_classic()+
        xlab('overlap proportion mean')+
        ggtitle(paste0('chromosome Y, p=', round(mean(v2_pm >v2),3)))
      pdf(paste0(pdir, path, '/gender/gender_hist_mean_of_overlap_proportion_PmReportedGeneOrder_', i, '.pdf'), width = 5, height = 2.6)
      gridExtra::grid.arrange(p1,p2,nrow=1)
      dev.off()
    }
  }
}

# ### compare with limma
# res.lm <- readRDS(paste0('hca/real/testvar/result/limma/', path, '/gender_res.rds'))
# allg.lm <- sub(':.*', '', rownames(res.lm))
# v1.lm <- cumsum(allg.lm %in% u1)/seq(1,length(allg.lm))
# v2.lm <- cumsum(allg.lm %in% u2)/seq(1,length(allg.lm))  
# saveRDS(mean(v1.lm), paste0('hca/real/testvar/result/limma/', path, '/gender_chrX_overlap.rds'))
# saveRDS(mean(v2.lm), paste0('hca/real/testvar/result/limma/', path, '/gender_chrY_overlap.rds'))


