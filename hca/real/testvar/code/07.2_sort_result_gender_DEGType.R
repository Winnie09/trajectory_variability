library(here)
library(ggplot2)
setwd(here())
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'

for (path in c('erythroid', 'lymph', 'monocyte')){
  # ---------
  # show data
  # ---------
  ## read data
  Res <- readRDS(paste0(ddir, path, '/gender_res.rds'))
  res <- Res$statistics
  ## identify DEGType 
  DEGType <- getDEGType(Res)
  res <- cbind(res, DEGType = DEGType[rownames(res)])
  saveRDS(res, paste0(rdir, path, '/gender_fdr_res.rds'))
  write.csv(res, paste0(rdir, path, '/gender_fdr_res.csv'))
  sink(paste0(rdir, path, '/gender_DEGType.txt'))
  table(DEGType)
  sink()
  ## plot DEG: all, meanSig, trendSig
  res <- res[order(res[,7], -res[,8]), ]
  pdf(paste0(pdir, path, '/gender_diffgene.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, rownames(res)[1:min(25,nrow(res))], type = 'variable', sep = ':.*')
  dev.off()
  
  res <- res[order(res[,1], -res[,2]), ]
  int <- intersect(rownames(res)[res[,1]<0.05], names(DEGType)[DEGType == 'meanSig'])
  pdf(paste0(pdir, path, '/gender_diffgene_meanSig.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
  dev.off()
  
  res <- res[order(res[,4], -res[,5]), ]
  int <- intersect(rownames(res)[res[,4]<0.05], names(DEGType)[DEGType == 'trendSig'])
  pdf(paste0(pdir, path, '/gender_diffgene_trendSig.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
  dev.off()
  
  
  # ------------
  #  Evaluation 
  # ------------
  ## read in gold standard Sex difference genes (chrX, chrY)
  u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
  u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
  ## TP
  res <- res[order(res[,7], -res[,8]), ]
  allg <- sub(':.*', '', rownames(res))
  diffgene = sub(':.*', '', rownames(res[res[,7] < 0.05,]))
  str(diffgene)
  
  diffgene.full <- rownames(res[res[,7] < 0.05,])
  diffgene.full <- diffgene.full[diffgene %in% u1]
  pdf(paste0(pdir, path, '/gender_true_diffgene_chrX.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, diffgene.full[1:min(25,length(diffgene.full))], type = 'variable', sep = ':.*')
  dev.off()
  
  diffgene.full <- rownames(res[res[,7] < 0.05,])
  diffgene.full <- diffgene.full[diffgene %in% u2]
  pdf(paste0(pdir, path, '/gender_true_diffgene_chrY.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, diffgene.full[1:min(25,length(diffgene.full))], type = 'variable', sep = ':.*')
  dev.off()
  
  ### proportion with Sex gold standard
  #### using all DEGs
  v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
  v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
  saveRDS(mean(v1), paste0(rdir, path, '/gender_chrX_overlap.rds'))
  saveRDS(mean(v2), paste0(rdir, path, '/gender_chrY_overlap.rds'))
  
  v1_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(allg, length(intersect(allg,u1)))
    v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
  },mc.cores=detectCores()))
  summary(colMeans(v1_pm))
  
  v2_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(allg, length(intersect(allg,u2)))
    v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
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
    ggtitle('chrX')
  
  p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
   geom_density(alpha=.2, fill="lightblue") +
   geom_vline(xintercept = mean(v1), color = 'darkblue') +
    theme_classic()+
    xlab('overlap proportion mean')+
    ggtitle('chrY')
  
  pdf(paste0(pdir, path, '/gender_hist_overlap_proportion_mean.pdf'), width = 7, height = 3)
  gridExtra::grid.arrange(p1,p2,nrow=1)
  dev.off()
  
  ## add here: the overlap using only trendSig, meanSig, other, bothSig
  for (i in unique(DEGType)){
    print(i)
    trendSig <- names(DEGType)[DEGType == i]
    allg <- sub(':.*','',intersect(rownames(res[order(res[,4], -res[,5]), ]), trendSig))
    v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
    v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
    saveRDS(mean(v1), paste0(rdir, path, '/gender_chrX_overlap_',i,'.rds'))
    saveRDS(mean(v2), paste0(rdir, path, '/gender_chrY_overlap_',i,'.rds'))
    
    v1_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
      set.seed(myseed)
      w1 = sample(allg, length(intersect(allg,u1)))
      v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
    },mc.cores=detectCores()))
    summary(colMeans(v1_pm))
    
    v2_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
      set.seed(myseed)
      w1 = sample(allg, length(intersect(allg,u2)))
      v1 <- cumsum(allg %in% w1)/seq(1,length(allg))
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
      ggtitle('chrX')
    
    p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
     geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
     geom_density(alpha=.2, fill="lightblue") +
     geom_vline(xintercept = mean(v1), color = 'darkblue') +
      theme_classic()+
      xlab('overlap proportion mean')+
      ggtitle('chrY')
    
    pdf(paste0(pdir, path, '/gender_hist_overlap_proportion_mean',i,'.pdf'), width = 7, height = 3)
    print(gridExtra::grid.arrange(p1,p2,nrow=1))
    dev.off()
  
  }
  
  ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  ###
  df = data.frame(chrX=v1, chrY=v2, chrX_pm = rowMeans(v1_pm), chrY_pm = rowMeans(v2_pm), order = seq(1,length(v1)))
  # saveRDS(df, './hca/geneexpr/result/df_chrX_chrY_pm_order.rds')
  
  ###
  mat <- NULL
  for (i in 1:4) {
    mat <- rbind(mat,data.frame(v=df[,i],order=df[,5],type=colnames(df)[i]))
  }
  library(ggplot2)
  library(RColorBrewer)
  pdf(paste0(pdir, path, '/gender_curve_number.pdf'), width=3.5, height=3)
  print(ggplot(mat,aes(x=order,y=v,col=type, fill=type), alpha=.2) + 
    geom_line() + 
    xlim(c(0,30)) + 
    ylim(c(0,10))+
    theme_classic()+
    ylab('number of ChrX/Y genes') + 
    xlab('top n genes')+
    scale_color_manual(values = brewer.pal(14, 'Paired')[c(1,7,2,8)]))
  dev.off()
  
  ### compare with limma
  res.lm <- readRDS(paste0(ddir, path,'/meandiff_gender_res.rds'))
  allg.lm <- sub(':.*', '', rownames(res.lm))
  v1.lm <- cumsum(allg.lm %in% u1)/seq(1,length(allg.lm))
  v2.lm <- cumsum(allg.lm %in% u2)/seq(1,length(allg.lm))  
  saveRDS(mean(v1.lm), paste0(rdir, path, '/meandiff_gender_chrX_overlap.rds'))
  saveRDS(mean(v2.lm), paste0(rdir, path, '/meandiff_gender_chrY_overlap.rds'))
  
}




