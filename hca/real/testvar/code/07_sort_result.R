library(here)
setwd(here())
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'

for (path in c('erythroid', 'lymph', 'monocyte')){
  dir.create(paste0(pdir, path))
  Res <- readRDS(paste0(ddir, path, '/age_res.rds'))
  res = data.frame(fdr = Res$fdr, pvalue = Res$pvalue, fc = Res$foldchange,
                   stringsAsFactors = FALSE)
  rownames(res) <- names(Res$fdr)
  res = res[order(res[,1], -abs(res[,3])),]
  saveRDS(res, paste0(rdir, path, '/age_fdr_res.rds'))
  write.csv(res, paste0(rdir, path, '/age_fdr_res.csv'))
  
  res <- res[res[,1]<0.05,]
  pdf(paste0(pdir, path, '/age_diffgene.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, rownames(res)[1:min(25, nrow(res))], variable = 'age', sep = ':.*')
  dev.off()
  
  ###
  Res <- readRDS(paste0(ddir, path, '/gender_res.rds'))
  res = data.frame(fdr = Res$fdr, pvalue = Res$pvalue, fc = Res$foldchange,
                   stringsAsFactors = FALSE)
  rownames(res) <- names(Res$fdr)
  res = res[order(res[,1], -abs(res[,3])),]
  saveRDS(res, paste0(rdir, path, '/gender_fdr_res.rds'))
  write.csv(res, paste0(rdir, path, '/gender_fdr_res.csv'))
  
  pdf(paste0(pdir, path, '/gender_diffgene.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, rownames(res)[1:min(25,nrow(res))], variable = 'gender', sep = ':.*')
  dev.off()
  
  u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
  u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
  allg = sub(':.*', '', rownames(res[res[,1] < 0.05,]))
  str(allg)
  
  allg.full <- rownames(res[res[,1] < 0.05,])
  allg.full <- allg.full[allg %in% u1]
  pdf(paste0(pdir, path, '/gender_true_diffgene_chrX.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, allg.full[1:min(25,length(allg.full))], variable = 'gender', sep = ':.*')
  dev.off()
  
  allg.full <- rownames(res[res[,1] < 0.05,])
  allg.full <- allg.full[allg %in% u2]
  pdf(paste0(pdir, path, '/gender_true_diffgene_chrY.pdf'), width = 12, height = 12)
  plotGenePopulation(Res, allg.full[1:min(25,length(allg.full))], variable = 'gender', sep = ':.*')
  dev.off()

  v1 <- cumsum(allg %in% u1)
  v2 <- cumsum(allg %in% u2)
  
  v1_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(intersect(allg,u1), length(intersect(allg,u1)), replace = TRUE)
    v1 <- cumsum(allg %in% w1)
  },mc.cores=detectCores()))
  summary(colMeans(v1_pm))
  
  v2_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(intersect(allg,u2), length(intersect(allg,u2)), replace = TRUE)
    v1 <- cumsum(allg %in% w1)
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
    xlab('mean of cumulated sum')+
    ggtitle('chrX')
  
  p2 <- ggplot(data = pd, aes(x=chrY_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white", stat = )+
   geom_density(alpha=.2, fill="lightblue") +
   geom_vline(xintercept = mean(v1), color = 'darkblue') +
    theme_classic()+
    xlab('mean of cumulated sum')+
    ggtitle('chrY')
  pdf(paste0(pdir, path, '/gender_hist_mean_cumulated_sum.pdf'), width = 7, height = 3)
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
}

