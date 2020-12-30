library(here)
setwd(here())
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'
u1 <- readRDS('agediff/result/age_diffgene.rds')
source('function/01_function.R')

for (path in c('monocyte', 'lymph', 'erythroid')){
  print(path)
  Res <- readRDS(paste0(ddir, path, '/age_res.rds'))
  res <- readRDS(paste0('hca/real/testvar/result/', path, '/age_fdr_res.rds'))

  allg = sub(':.*', '', rownames(res[res[,1] < 0.05,]))
  str(allg)
  
  v1 <- cumsum(allg %in% u1)
  
  library(parallel)
  v1_pm <- do.call(cbind,mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    w1 = sample(intersect(allg,u1), length(intersect(allg,u1)), replace = TRUE)
    v1 <- cumsum(allg %in% w1)
  },mc.cores=detectCores()))
  summary(colMeans(v1_pm))
  rownames(v1_pm) <- paste0('top',seq(1,nrow(v1_pm)))
  
  pd <- data.frame(agediff_pm = colMeans(v1_pm), stringsAsFactors = FALSE)
  library(ggplot2)
  pdf(paste0(pdir, path, '/agediff_hist_mean_cumulated_sum.pdf'), width = 3.5, height = 3)
  print(ggplot(data = pd, aes(x=agediff_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white")+
   geom_density(alpha=.2, fill="#FF6666") +
   geom_vline(xintercept = mean(v1), color = 'red') +
    theme_classic()+
    xlab('mean of cumulated sum')+
    ggtitle('age differential'))
  dev.off()
  
  df = data.frame(ourmethod=v1, agediff_pm = rowMeans(v1_pm), order = seq(1,length(v1)))
  
  mat <- NULL
  for (i in 1:2) {
    mat <- rbind(mat,data.frame(v=df[,i],order=df[,3],type=colnames(df)[i]))
  }
  library(ggplot2)
  library(RColorBrewer)
  pdf(paste0(pdir, path, '/agediff_curve_number.pdf'), width=3.5, height=3)
  print(ggplot(mat,aes(x=order,y=v,col=type, fill=type), alpha=.2) + 
    geom_line() + 
    xlim(c(0,30)) + 
    ylim(c(0,10))+
    theme_classic()+
    ylab('#GTEX age genes') + 
    xlab('top n genes')+
    scale_color_manual(values = brewer.pal(14, 'Paired')[c(2,8)]))
  dev.off()
  
  pdf(paste0(pdir, path, '/agediff_true_genes.pdf'), width=12, height=12)
  allg.full <- rownames(res[res[,1] < 0.05,])
  allg.full <- allg.full[allg %in% u1]
  plotGenePopulation(Res, allg.full[1:min(length(allg.full), 25)], variable = 'age')
  dev.off()
  
}



