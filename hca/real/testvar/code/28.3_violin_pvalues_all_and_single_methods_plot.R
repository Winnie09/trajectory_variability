rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
library(parallel)
# setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home/whou10/scratch16/whou10/')
# source('function/01_function.R')
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)

## "erythroid" "lymph"     "monocyte"
for (celltype in setdiff(list.files('trajectory_variability/hca/real/testvar/result/EM_pm'),'perf')) {
  res = readRDS(paste0('trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin_plotdata.rds'))
  permud = res[[1]]
  reald = res[[2]]
  exclude = c('Lamian_excTradeSeqPT', 'Lamian_excTradeSeqET', 'Lamian_excTradeSeqDT', 'tradeSeqET', 'tradeSeqPT', 'tradeSeqDT', 'tradeSeqET_excludeLamian', 'tradeSeqPT_excludeLamian', 'tradeSeqDT_excludeLamian', 'monocle2trajTest', 'monocle2trajTest_excLamian', 'Lamian_excMonocle2trajTest')
  permud = permud[!permud[,3] %in% exclude, ]
  reald = reald[!reald[,4] %in% exclude, ]
  
  tmp = reald[reald[,2] == 'chrX',]
  xv = tmp[,1]
  names(xv) = tmp[,4] 
  allm = names(sort(xv))
  
  m1 = allm[grepl('excLamian', allm)]
  m2 = allm[grepl('Lamian_exc', allm)]
  m3 =  setdiff(allm, c(m1, m2))
  permud[,3] = factor(permud[,3], levels = c(m1, m2, m3))
  reald[,4] = factor(reald[,4], levels = c(m1, m2, m3))
  
  
  pdf(paste0('trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin_all.pdf'),width=11,height=3)
  print(ggplot() + 
          geom_violin(data=permud,aes(x=method,y=per,col=type)) + 
          geom_point(data=reald,aes(x=method,y=per,col=type),size=1) + 
          geom_text(data=reald,aes(x=method,y=max(reald$per)*1.3,label=pvalue), size = 10*5/14)+
          #theme_compact() + 
          facet_wrap(~type) + 
          coord_flip(ylim=c(0,max(reald$per)*1.5)) + 
          xlab('') + 
          ylab('Proportion') + 
          scale_color_manual(values=c('chrX'=brewer.pal(3,'Pastel1')[1],'chrY'=brewer.pal(3,'Pastel1')[2])) +
          theme(legend.position = 'bottom',strip.background = element_blank(),strip.text = element_text(size=7),legend.title = element_blank(), text = element_text(size = 5)) + 
          scale_y_continuous(breaks=c(0,round(max(reald$per)*0.6,2),round(max(reald$per)*1.2,2)))
  )
  dev.off()
  
  pdf(paste0('trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin_single_method.pdf'),width=6,height=1.5)
  print(ggplot() + 
          geom_violin(data=permud[permud[,3]%in%m3, ],aes(x=method,y=per,col=type)) + 
          geom_point(data=reald[reald[,4]%in%m3, ],aes(x=method,y=per,col=type),size=1) + 
          geom_text(data=reald[reald[,4]%in%m3, ],aes(x=method,y=max(reald$per)*1.3,label=pvalue), size = 5)+
          #theme_compact() + 
          facet_wrap(~type) + 
          coord_flip(ylim=c(0,max(reald$per)*1.5)) + 
          xlab('') + 
          ylab('Proportion') + 
          scale_color_manual(values=c('chrX'=brewer.pal(3,'Pastel1')[1],'chrY'=brewer.pal(3,'Pastel1')[2])) +
          theme(legend.position = 'bottom',strip.background = element_blank(),strip.text = element_text(size=7), legend.title = element_blank(), text = element_text(size = 5)) + 
          scale_y_continuous(breaks=c(0,round(max(reald$per)*0.6,2),round(max(reald$per)*1.2,2)))
  )
  dev.off()
}
