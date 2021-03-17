library(here)
library(ggplot2)
setwd(here())
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/EM_pm/'
pdir <- 'hca/real/testvar/plot/EM_pm/'

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')


for (path in c('erythroid', 'lymph', 'monocyte')){
  ## read data
  Res <- readRDS(paste0(ddir, path, '/gender/gender_res.rds'))
  res <- Res$statistics
  ## identify DEGType 
  DEGType <- getDEGType(Res)
  # ------------------
  #  Evaluation Gender 
  # ------------------
  ## TP
  res <- res[order(res[,'fdr.overall'], -res[,'z.overall']), ]
  allg <- sub(':.*', '', rownames(res))
  diffgene = sub(':.*', '', rownames(res[res[,'fdr.overall'] < 0.05,]))
  str(diffgene)
  if (length(diffgene) > 0) {
    diffgene.full <- rownames(res[res[,'fdr.overall'] < 0.05,])
    diffgene.full <- diffgene.full[diffgene %in% u1]
    if (!is.na(diffgene.full)){
      pdf(paste0(pdir, path, '/gender/gender_true_diffgene_chrX.pdf'), width = 12, height = 12)
      plotGenePopulation(Res, diffgene.full[1:min(25,length(diffgene.full))], type = 'variable', sep = ':.*')
      dev.off()
    }
    diffgene.full <- diffgene.full[diffgene %in% u2]
    if (!is.na(diffgene.full)){
      pdf(paste0(pdir, path, '/gender_true_diffgene_chrY.pdf'), width = 12, height = 12)
      plotGenePopulation(Res, diffgene.full[1:min(25,length(diffgene.full))], type = 'variable', sep = ':.*')
      dev.off()
    }
  }
  ### proportion with Sex gold standard
  #### using all DEGs
  v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
  v1 <- mean(v1[seq(1,length(v1)) %% 10 == 0])
  v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
  v2 <- mean(v2[seq(1,length(v2)) %% 10 == 0])
  saveRDS(v1, paste0(rdir, path, '/gender/gender_chrX_overlap.rds'))
  saveRDS(v2, paste0(rdir, path, '/gender/gender_chrY_overlap.rds'))
  
  ##### permute gold standard gene sets
  v1_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(allg, length(intersect(allg,u1)))
    v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
    mean(v1[seq(1,length(v1)) %% 10 == 0])
  },mc.cores=detectCores()))
  
  
  v2_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(allg, length(intersect(allg,u2)))
    v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
    mean(v1[seq(1,length(v1)) %% 10 == 0])
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
    ggtitle(paste0('chrX, p=', round(mean(v1_pm > v1),3)))
  
  p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
    geom_density(alpha=.2, fill="lightblue") +
    geom_vline(xintercept = v2, color = 'darkblue') +
    theme_classic()+
    xlab('overlap proportion mean')+
    ggtitle(paste0('chrY, p=', round(mean(v2_pm >v2),3)))
  
  pdf(paste0(pdir, path, '/gender/gender_hist_mean_of_overlap_proportion.pdf'), width = 7, height = 3)
  gridExtra::grid.arrange(p1,p2,nrow=1)
  dev.off()
  
  #### using only one type of DEG: trendSig, meanSig, other, bothSig
  if (length(diffgene) > 0){
    for (i in setdiff(unique(DEGType), 'nonDEG')){
      print(i)
      SigGene <- names(DEGType)[DEGType == i]
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
      saveRDS(v1, paste0(rdir, path, '/gender/gender_chrX_overlap_', i,'.rds'))
      saveRDS(v2, paste0(rdir, path, '/gender/gender_chrY_overlap_', i, '.rds'))
      
      ##### permute gold standard genes
      v1_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
        set.seed(myseed)
        w1 = sample(allg, length(intersect(allg,u1)))
        v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
        mean(v1[seq(1,length(v1)) %% 10 == 0])
      },mc.cores=detectCores()))
      
      
      v2_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
        set.seed(myseed)
        w1 = sample(allg, length(intersect(allg,u2)))
        v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
        mean(v1[seq(1,length(v1)) %% 10 == 0])
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
        ggtitle(paste0('chrX, p=', round(mean(v1_pm > v1),3)))
      
      p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
        geom_density(alpha=.2, fill="lightblue") +
        geom_vline(xintercept = v2, color = 'darkblue') +
        theme_classic()+
        xlab('overlap proportion mean')+
        ggtitle(paste0('chrY, p=', round(mean(v2_pm >v2),3)))
      pdf(paste0(pdir, path, '/gender/gender_hist_mean_of_overlap_proportion_', i, '.pdf'), width = 7, height = 3)
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
