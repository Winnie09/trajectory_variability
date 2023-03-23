rm(list=ls())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
library(here)
setwd(here())
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/perf/'
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr_pm_window/plot/perf/'
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/data/selgene/selgene1.rds') # trend diff only
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/data/selgene/selgene2.rds') # mean diff only
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/data/data/selgene/selgene3.rds') # trend and mean
selgene = selgene1
rmgene = c(selgene2, selgene3)
str(selgene)
str(rmgene)
dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
source('function/01_function.R')
library(gridExtra)

perflist <- list()
m = 'Lamian.pm'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.trendDiff'], zscore = r$statistics[,'z.trendDiff'], stringsAsFactors = F)
  res = res[complete.cases(res),]
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'trend.test', AreaUnderSensFdr(SensFdr(c(selgene1,selgene3), res)), "Lamian")
})
df1.nc[3,2] = as.numeric(df1.nc[3,1])*0.455+as.numeric(df1.nc[3,3])*0.558
df1.nc[4,2] = as.numeric(df1.nc[4,1])*0.498+as.numeric(df1.nc[4,3])*0.502
perflist[['Lamian.trend.test']] = t(df1.nc)

df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.meanDiff'], zscore = r$statistics[,'z.meanDiff'], stringsAsFactors = F)
  res = res[complete.cases(res),]
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'mean.test', AreaUnderSensFdr(SensFdr(c(selgene2,selgene3), res)), 'Lamian')
})
df1.nc[3,4] = (as.numeric(df1.nc[3,3])+as.numeric(df1.nc[3,5]))/2
df1.nc[4,4] = (as.numeric(df1.nc[4,3])+as.numeric(df1.nc[4,5]))/2
perflist[['Lamian.mean.test']] = t(df1.nc)

###
af = list.files(paste0(ddir, 'Lamian.chisq/'))
af <- af[grepl('.rds', af)]
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, 'Lamian.chisq/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.chisq.trendDiff'], stat = (r$ll3-r$ll2)/r$statistics[,'df.diff.trendDiff'], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'trend.test', AreaUnderSensFdr(SensFdr(c(selgene1, selgene3), res)), 'Lamian.chisq')
})
perflist[['Lamian.chisq.trend.test']] = t(df1)

df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, 'Lamian.chisq/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.chisq.meanDiff'], stat = (r$ll2-r$ll1)/r$statistics[,'df.diff.meanDiff'], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'mean.test', AreaUnderSensFdr(SensFdr(c(selgene2, selgene3), res)), 'Lamian.chisq')
})
perflist[['Lamian.chisq.mean.test']] = t(df1)


###
perf <- do.call(rbind,perflist)
perf <- data.frame('SignalStrength'=as.numeric(perf[,1]),'Comparison'=perf[,2],'Fdr.Diff'=as.numeric(perf[,3]),'AUC'=as.numeric(perf[,4]), Method = perf[,5])
saveRDS(perf, paste0(rdir,'perf_trenddiff_meandiff.rds'))


pdonly <- perf
# pdonly$comparison <- pdonly$Method
# pdonly$Method <- 'Lamian'

pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/perf/perf.rds')
pd[,2] <- as.character(pd[,2])
pd[pd[,2]=='Lamian_chisq',2] <- 'Lamian.chisq'
pd0 <- pd
pd0$Comparison <- 'overall.test'


pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/perf/perf_meanOnly.rds')
pd[,2] <- as.character(pd[,2])
pd1 <- pd
pd1$Comparison = 'mean only'


pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/perf/perf_trendOnly.rds')
pd[,2] <- as.character(pd[,2])
pd[pd[,2]=='Lamian_chisq',2] <- 'Lamian.chisq'
pd2 <- pd
pd2$Comparison = 'trend only'


pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr_pm_window/result/perf/perf_trendMean.rds')
pd[,2] <- as.character(pd[,2])
pd3 <- pd
pd3$Comparison = 'trend & mean'

pd <- rbind(pd0, pd1,pd2,pd3)
pd$Comparison <- factor(pd$Comparison, levels = c('overall.test', 'trend only', 'trend & mean','mean only'))


pd <- rbind(pd,pdonly)
pd <- pd[complete.cases(pd), ]
library(ggplot2)
library(RColorBrewer)
colv = c(brewer.pal(9,'Set1')[c(1:2)], brewer.pal(8, 'Dark2')[1:8])
names(colv) <- c('Lamian', 'Lamian.chisq', setdiff(unique(pd[,2]), c('Lamian', 'Lamian.chisq')))


# pdf(paste0(pdir, 'compare_fdr_diff_auc_trenddiff_meandiff.pdf'),width=12.5,height=2.8) ## good 
pdf(paste0(pdir, 'compare_fdr_diff_auc_trenddiff_meandiff.pdf'),width=16,height=4) ##  large
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=1)  + 
  geom_line(size=0.5) + 
  theme_bw() + 
  xlim(c(0,4))+
  ylab('FDR.difference') +
scale_color_manual(values = colv)  + 
  facet_wrap(~Comparison, nrow =1) +
  theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(-0.04,0.2))+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=1)  + 
  geom_line(size=0.5) + 
  theme_bw() + 
    xlim(c(0,4))+
    ylab('AUC (sensitivity-realFDR)') +
  scale_color_manual(values = colv)  + 
  facet_wrap(~Comparison, nrow =1)+
    theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(0,1))+
    guides(color=guide_legend(nrow=2,byrow=TRUE))
grid.arrange(p1,p2,nrow=1)
dev.off()

pd[,5] = factor(as.character(pd[,5]), levels = c("overall.test","trend.test", "mean.test", "trend only","trend & mean", "mean only"))
pdf(paste0(pdir, 'compare_fdr_diff_auc_trenddiff_meandiff_2row.pdf'),width=7,height=4.3)
# pdf(paste0(pdir, 'compare_fdr_diff_auc_trenddiff_meandiff_2row.pdf'),width=15,height=6.5)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=0.5)  + 
  geom_line(size=0.3) + 
  xlim(c(0,4))+
  ylab('FDR.difference') +
  scale_color_manual(values = colv)  + 
  theme_bw() +
  theme(strip.text.x = element_text(size = 10), axis.title = element_text(size = 10), axis.text = element_text(size = 10)) + 
  facet_wrap(~Comparison, nrow =2) +
  theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(-0.04,0.2))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, label.position = "right",
                             title.position = "left", title.vjust = 1, legend.spacing.x = unit(0.1, 'cm'),
                            legend.spacing.y = unit(0.1, 'pt'), legend.key.size = unit(0.5, "cm")))
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=0.5)  + 
  geom_line(size=0.3) + 
  xlim(c(0,4))+
  ylab('AUC (sensitivity-realFDR)') +
  scale_color_manual(values = colv)  + 
  theme_bw() +
  theme(strip.text.x = element_text(size = 10), axis.title = element_text(size = 10), axis.text = element_text(size = 10)) + 
  facet_wrap(~Comparison, nrow =2)+
  theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(0,1))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, label.position = "right",
                             title.position = "left", title.vjust = 1, legend.spacing.x = unit(0.1, 'cm'),
                            legend.spacing.y = unit(0.1, 'pt')))
grid.arrange(p1,p2,nrow=1)
dev.off()


