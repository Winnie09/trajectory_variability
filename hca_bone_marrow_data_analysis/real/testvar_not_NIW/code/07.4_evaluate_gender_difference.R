rm(list=ls())
library(here)
library(ggplot2)
setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# setwd(here())
source('function/01_function.R')
u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')
  
  
for (path in c('erythroid', 'lymph', 'monocyte')){
  ddir <- rdir <- paste0('hca/real/testvar/result/', path)
  pdir <- paste0('hca/real/testvar/plot/', path)
  dir.create(pdir)
  ### ------------------------------
  ### gender difference: Evaluation 
  ### ------------------------------
  Res <- readRDS(paste0(ddir, '/gender/gender_res.rds'))
  res <- Res$statistics
  DDGType <- getDDGType(Res)
  res <- cbind(res, DDGType = DDGType[rownames(res)])
  ## read in gold standard Sex difference genes (chrX, chrY)    
  res <- res[order(res[,7], -res[,8]), ]
  allg <- sub(':.*', '', rownames(res))
  diffgene = sub(':.*', '', rownames(res[res[,7] < 0.05,]))
  str(diffgene)
  
  diffgene.full <- rownames(res[res[,7] < 0.05,])
  diffgene.full <- diffgene.full[diffgene %in% u1]
  if (length(diffgene.full) > 0){
    png(paste0(pdir, '/gender/gender_true_diffgene_chrX.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, diffgene.full[1:min(25,length(diffgene.full))], type = 'variable', sep = ':.*')
    dev.off()
  }
  
  diffgene.full <- rownames(res[res[,7] < 0.05,])
  diffgene.full <- diffgene.full[diffgene %in% u2]
  if (length(diffgene.full) > 0){
    png(paste0(pdir, '/gender/gender_true_diffgene_chrY.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, diffgene.full[1:min(25,length(diffgene.full))], type = 'variable', sep = ':.*')
    dev.off()
  }
  
  ### overlap proportion with Sex gold standard
  v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
  v1 <- v1[seq(1,length(v1)) %% 10 == 0]
  v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
  v2 <- v2[seq(1,length(v2)) %% 10 == 0]
  saveRDS(mean(v1), paste0(rdir, path, 'gender/gender_chrX_overlap.rds'))
  saveRDS(mean(v2), paste0(rdir, path, 'gender/gender_chrY_overlap.rds'))
  
  v1_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed+100)
    allg.pm = sample(allg)
    tmp = cumsum(allg.pm %in% u1)/seq(1,length(allg.pm))
    tmp <- tmp[seq(1,length(tmp)) %% 10 == 0]
  },mc.cores=detectCores()))
  summary(colMeans(v1_pm))
  mean(v1_pm>mean(v1))
  
  v2_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    allg.pm = sample(allg)
    tmp = cumsum(allg.pm %in% u2)/seq(1,length(allg.pm))
    tmp <- tmp[seq(1,length(tmp)) %% 10 == 0]
  },mc.cores=detectCores()))
  summary(colMeans(v2_pm))
  
  rownames(v1_pm) <- paste0('top',seq(1,nrow(v1_pm)))
  rownames(v2_pm) <- paste0('top',seq(1,nrow(v2_pm)))
  
  pd <- data.frame(chrX_pm = colMeans(v1_pm), chrY_pm = colMeans(v2_pm), stringsAsFactors = FALSE)
  p1 <- ggplot(data = pd, aes(x=chrX_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white")+
   geom_density(alpha=.2, fill="#FF6666") +
   geom_vline(xintercept = mean(v1), color = 'red') +
    theme_classic()+
    xlab('overlap proportion mean')+
    ggtitle(paste0('chrX, p.value=', round(mean(v1_pm>mean(v1)),3)))
  
  p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
   geom_density(alpha=.2, fill="lightblue") +
   geom_vline(xintercept = mean(v1), color = 'darkblue') +
    theme_classic()+
    xlab('overlap proportion mean')+
    ggtitle(paste0('chrY, p.value=', round(mean(v2_pm>mean(v2)),3)))
  
  pdf(paste0(pdir, '/gender/gender_hist_overlap_proportion_mean.pdf'), width = 6, height = 2.6)
  gridExtra::grid.arrange(p1,p2,nrow=1)
  dev.off()
  
  df = data.frame(chrX=v1, chrY=v2, chrX_pm = rowMeans(v1_pm), chrY_pm = rowMeans(v2_pm), order = seq(1,length(v1)))
  # saveRDS(df, './hca/geneexpr/result/df_chrX_chrY_pm_order.rds')
  
  mat <- NULL
  for (i in 1:4) {
    mat <- rbind(mat,data.frame(v=df[,i],order=df[,5],type=colnames(df)[i]))
  }
  library(ggplot2)
  library(RColorBrewer)
  pdf(paste0(pdir, '/gender/gender_curve_number.pdf'), width=3.5, height=3)
  ggplot(mat,aes(x=order,y=v,col=type, fill=type), alpha=.2) + 
    geom_line() + 
    xlim(c(0,30)) + 
    ylim(c(0,max(mat[,1])+0.2))+
    theme_classic()+
    ylab('number of ChrX/Y genes') + 
    xlab('top n genes')+
    scale_color_manual(values = brewer.pal(14, 'Paired')[c(2,1,8,7)])
  dev.off()
  
  ### compare with limma
  res.lm <- readRDS(paste0(ddir,'/gender/meandiff_gender_res.rds'))
  allg.lm <- sub(':.*', '', rownames(res.lm))
  v1.lm <- cumsum(allg.lm %in% u1)/seq(1,length(allg.lm))
  v1.lm <- v1.lm[seq(1,length(v1.lm)) %% 10 == 0]
  v2.lm <- cumsum(allg.lm %in% u2)/seq(1,length(allg.lm))  
  v2.lm <- v2.lm[seq(1,length(v2.lm)) %% 10 == 0]
  saveRDS(mean(v1.lm), paste0(rdir, '/gender/meandiff_gender_chrX_overlap.rds'))
  saveRDS(mean(v2.lm), paste0(rdir, '/gender/meandiff_gender_chrY_overlap.rds'))
}

