library(here)
setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# setwd(here())
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'
u1 <- readRDS('agediff/result/age_diffgene.rds')
source('function/01_function.R')
res.gtex <- readRDS('agediff/result/res.rds')

for (path in c('monocyte', 'lymph', 'erythroid')){
  print(path)
  res <- readRDS(paste0('hca/real/testvar/result/', path, '/age/age_fdr_res.rds'))
  res <- res[order(res[,7], -abs(res[,8])), ]
  rn <- sub(':.*', '', rownames(res))
  res <- res[!duplicated(rn), ]
  rownames(res) <- rn[!duplicated(rn)]
  res <- res[rownames(res) %in% rownames(res.gtex), ]
  
  res.lm <- readRDS(paste0('hca/real/testvar/result/', path, '/age/meandiff_age_res.rds'))
  res.lm <- res.lm[order(res.lm[,5], -abs(res.lm[,1])), ]
  rn <- sub(':.*', '', rownames(res.lm))
  res.lm <- res.lm[!duplicated(rn), ]
  rownames(res.lm) <- rn[!duplicated(rn)]
  res.lm <- res.lm[rownames(res.lm) %in% rownames(res.gtex), ]
  
  ## our method overlap
  allg <- rownames(res)
  
  targetgene <- allg
  v1 <- cumsum(targetgene %in% intersect(targetgene,u1))/seq(1:length(targetgene))
  v1 <- v1[seq(1,length(v1)) %% 10 == 0]
  saveRDS(mean(v1), paste0(rdir, path, '/age/age_overlap.rds'))
  
  ## limmma overlap
  targetgene2 <- rownames(res.lm)
  v2 <- cumsum(targetgene2 %in% intersect(targetgene2,u1))/seq(1:length(targetgene2))
  v2 <- v2[seq(1,length(v2)) %% 10 == 0]
  saveRDS(mean(v2), paste0(rdir, path, '/age/meandiff_age_overlap.rds'))
  
  v1_pm <- sapply(seq(1,1e4), function(myseed){
    set.seed(myseed)
    sampgene <- sample(targetgene)
    k1 <- cumsum(sampgene %in% intersect(sampgene,u1))/seq(1:length(sampgene))
    k1 <- k1[seq(1,length(k1)) %% 10 == 0]
  })
  summary(colMeans(v1_pm))
  rownames(v1_pm) <- paste0('top',seq(1,nrow(v1_pm)))
  pval <- mean(v1_pm > mean(v1))
  
  pd <- data.frame(agediff_pm = colMeans(v1_pm), stringsAsFactors = FALSE)
  library(ggplot2)
  pdf(paste0(pdir, path, '/age/agediff_hist_mean_overlap_proportion.pdf'), width = 2.8, height = 2.5)
  print(ggplot(data = pd, aes(x=agediff_pm)) + 
   geom_histogram(aes(y=..density..), colour="black", fill="white")+
   geom_density(alpha=.2, fill="#FF6666") +
   geom_vline(xintercept = mean(v1), color = 'red') +
    theme_classic()+
    xlab('overlap proportion mean')+
    ggtitle(paste0('age, p.value = ', round(pval,3))))
  dev.off()
  
  
  df = data.frame(ourmethod=v1, agediff_pm = rowMeans(v1_pm), order = seq(1,length(v1)))
  mat <- NULL
  for (i in 1:2) {
    mat <- rbind(mat,data.frame(v=df[,i],order=df[,3],type=colnames(df)[i]))
  }
  library(ggplot2)
  library(RColorBrewer)
  pdf(paste0(pdir, path, '/age/agediff_curve_number.pdf'), width=3.5, height=2.2)
  print(ggplot(mat,aes(x=order,y=v,col=type, fill=type), alpha=.2) + 
    geom_line() + 
    theme_classic()+
    ylab('#GTEX age genes') + 
    xlab('top n genes')+
    scale_color_manual(values = brewer.pal(14, 'Paired')[c(2,8)]))
  dev.off()
  
  
  diffgene <- rownames(res[res[,7] < 0.05, ])
  str(diffgene)
  if (sum(diffgene%in% u1) > 0){
    Res <- readRDS(paste0(ddir, path, '/age/age_res.rds'))
    pdf(paste0(pdir, path, '/age/agediff_true_genes.pdf'), width=12, height=12)
    res <- readRDS(paste0('hca/real/testvar/result/', path, '/age/age_fdr_res.rds'))
    res <- res[order(res[,7], -abs(res[,8])), ]
    int <- intersect(rownames(res[res[,11] == 'TRUE', ]), rownames(res[res[,7] < 0.05, ]))
    plotGenePopulation(Res, int[1:min(length(int), 25)], type = 'variable')
    dev.off()
  }
}


