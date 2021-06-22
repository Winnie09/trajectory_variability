rm(list=ls())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
library(here)
setwd(here())
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds') # trend diff only
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds') # mean diff only
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds') # trend and mean
selgene = selgene1
rmgene = c(selgene2, selgene3)
str(selgene)
str(rmgene)
source('function/01_function.R')
library(gridExtra)

perflist <- list()
m = 'EM_pm'
af = list.files(paste0(ddir, m, '/'))
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.trendDiff'], zscore = r$statistics[,'z.trendDiff'], stringsAsFactors = F)
  res = res[complete.cases(res),]
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'trendDiff', AreaUnderSensFdr(SensFdr(c(selgene1,selgene3), res)))
})
perflist[['trendDiff']] = t(df1.nc)

df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(adj.P.Val = r$statistics[,'fdr.meanDiff'], zscore = r$statistics[,'z.meanDiff'], stringsAsFactors = F)
  res = res[complete.cases(res),]
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'meanDiff', AreaUnderSensFdr(SensFdr(c(selgene2,selgene3), res)))
})
# df1.nc[3,4] = (as.numeric(df1.nc[3,3]) + as.numeric(df1.nc[3,5]))/2
# df1.nc[4,4] = (as.numeric(df1.nc[4,3]) + as.numeric(df1.nc[4,5]))/2
perflist[['meanDiff']] = t(df1.nc)

perf <- do.call(rbind,perflist)
perf <- data.frame('SignalStrength'=as.numeric(perf[,1]),'Method'=perf[,2],'Fdr.Diff'=as.numeric(perf[,3]),'AUC'=as.numeric(perf[,4]))
saveRDS(perf, paste0(rdir,'perf_trenddiff_meandiff.rds'))


pdonly <- perf
pdonly$comparison <- sub('meanDiff','FDR.meanSig',sub('trendDiff','FDR.trendSig',pdonly$Method))
pdonly$Method <- 'Lamian'

pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf.rds')
pd[,2] <- as.character(pd[,2])
pd0 <- pd[pd[,2]!='EM_chisq', ]

pd0$comparison <- 'overall'
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_meanOnly.rds')
pd[,2] <- as.character(pd[,2])
pd1 <- pd[pd[,2]!='EM_chisq', ]

pd1$comparison = 'mean only'
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendOnly.rds')
pd[,2] <- as.character(pd[,2])
pd2 <- pd[pd[,2]!='EM_chisq', ]

pd2$comparison = 'trend only'
pd <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/result/perf/perf_trendMean.rds')
pd[,2] <- as.character(pd[,2])
pd3 <- pd[pd[,2]!='EM_chisq', ]

pd3$comparison = 'trend & mean'
pd <- rbind(pd0, pd1,pd2,pd3)
pd <- pd[pd[,2]!='TSCAN', ]
pd$comparison <- factor(pd$comparison, levels = c('overall', 'trend only', 'trend & mean','mean only'))
pd[pd[,2] == 'EM_pm', 2] <- 'Lamian' ## method name

pd <- rbind(pd,pdonly)
pd <- pd[complete.cases(pd), ]
library(ggplot2)
pdf(paste0(pdir, 'compare_fdr_diff_auc_trenddiff_meandiff.pdf'),width=12.5,height=2.8)
p1 <- ggplot(pd, aes(x = SignalStrength, y = Fdr.Diff, color=Method)) + 
  geom_point(size=2)  + 
  geom_line(size=1) + 
  theme_bw() + 
  xlim(c(0,4))+
  ylab('FDR.difference') +
  scale_color_brewer(palette = 'Set1')  +
  facet_wrap(~comparison, nrow =1) +
  theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(-0.04,0.2))+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p2 <- ggplot(pd, aes(x = SignalStrength, y = AUC, color=Method)) + 
  geom_point(size=2)  + 
  geom_line(size=1) + 
  theme_bw() + 
    xlim(c(0,4))+
    ylab('AUC (sensitivity-realFDR)') +
  scale_color_brewer(palette = 'Set1')  +
  facet_wrap(~comparison, nrow =1)+
    theme(legend.position = 'bottom',legend.title = element_blank()) + 
  ylim(c(0,1))+
    guides(color=guide_legend(nrow=2,byrow=TRUE))
grid.arrange(p1,p2,nrow=1)
dev.off()






