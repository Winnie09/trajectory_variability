rm(list=ls())
library(here)
setwd(here())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
selgene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene.rds')
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds')
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds')
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds')

source('function/01_function.R')
slist <- list()
m = 'EM_pm'
f = '4.rds'
print(f)
r = readRDS(paste0(ddir, m, '/', f))
res = data.frame(adj.P.Val = r$statistics[,'fdr.overall'], zscore = r$statistics[,'z.overall'], stringsAsFactors = F)
res = res[order(res[,1], res[,2]), ]
slist[['Lamian']] = rbind(cbind(SensFdr(selgene, res), method = 'Lamian', comparison = 'overall'),
cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'Lamian', comparison = 'trendOnly'),
cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'Lamian', comparison = 'meanOnly'),
cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'Lamian', comparison = 'trendMean'))

res = data.frame(adj.P.Val = r$statistics[,'fdr.trendDiff'], zscore = r$statistics[,'z.trendDiff'], stringsAsFactors = F)
res = res[order(res[,1], res[,2]), ]
res = res[complete.cases(res), ]
slist[['Lamian']] = rbind(slist[['Lamian']],cbind(SensFdr(c(selgene1,selgene3), res), method = 'Lamian', comparison = 'FDR.trendSig'))

res = data.frame(adj.P.Val = r$statistics[,'fdr.meanDiff'], zscore = r$statistics[,'z.meanDiff'], stringsAsFactors = F)
res = res[order(res[,1], res[,2]), ]
res = res[complete.cases(res), ]
slist[['Lamian']] = rbind(slist[['Lamian']],cbind(SensFdr(c(selgene1,selgene3), res), method = 'Lamian', comparison = 'FDR.meanSig'))

m = 'tradeSeq'
res = readRDS(paste0(ddir, m, '/', f))[[1]]
res = res[order(res[, 3],-abs(res[, 1])),]
slist[['tradeSeq_diffEndTest']] = rbind(cbind(SensFdr(selgene, res), method = 'tradeSeq_diffEndTest', comparison = 'overall'),
cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'tradeSeq_diffEndTest', comparison = 'trendOnly'),
cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'tradeSeq_diffEndTest', comparison = 'meanOnly'),
cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'tradeSeq_diffEndTest', comparison = 'trendMean'))

res = readRDS(paste0(ddir, m, '/', f))[[2]]
res[is.na(res[,3]),3] <- 1
res[is.na(res[,2]),2] <- 1
res[is.na(res[,1]),1] <- 0
res = res[order(res[, 3],-abs(res[, 1])),]
slist[['tradeSeq_patternTest']] = rbind(cbind(SensFdr(selgene, res), method = 'tradeSeq_patternTest', comparison = 'overall'),
cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'tradeSeq_patternTest', comparison = 'trendOnly'),
cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'tradeSeq_patternTest', comparison = 'meanOnly'),
cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'tradeSeq_patternTest', comparison = 'trendMean'))

res = readRDS(paste0(ddir, m, '/', f))[[3]]
res[is.na(res[,3]),3] <- 1
res[is.na(res[,2]),2] <- 1
res[is.na(res[,1]),1] <- 0
res = res[order(res[, 3],-abs(res[, 1])),]
slist[['tradeSeq_earlyDETest']] = rbind(cbind(SensFdr(selgene, res), method = 'tradeSeq_earlyDETest', comparison = 'overall'),
cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'tradeSeq_earlyDETest', comparison = 'trendOnly'),
cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'tradeSeq_earlyDETest', comparison = 'meanOnly'),
cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'tradeSeq_earlyDETest', comparison = 'trendMean'))

m = 'meandiff'
res = readRDS(paste0(ddir, m, '/', f))
res = res[order(res[,5], -abs(res[,1])), ]
slist[['limma']] = rbind(cbind(SensFdr(selgene, res), method = 'limma', comparison = 'overall'),
cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'limma', comparison = 'trendOnly'),
cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'limma', comparison = 'meanOnly'),
cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'limma', comparison = 'trendMean'))

## concatenate
saveRDS(slist, paste0(rdir,'sensfdr_list.rds'))
s <- do.call(rbind, slist)
str(s)
s <- as.data.frame(s)
s[,1] <- as.numeric(as.character(s[,1]))
s[,3] <- as.numeric(as.character(s[,3]))
s[,2] <- as.numeric(as.character(s[,2]))
s[,5] <- factor(as.character(s[,5]), levels = c("overall", "FDR.trendSig","FDR.meanSig","trendOnly","trendMean", "meanOnly" ))
saveRDS(s, paste0(rdir,'sensfdr.rds'))


##############################################
library(ggplot2)
library(ggthemes)
pdf(paste0(pdir, 'realFDR_reportedFDR_strength4.pdf'),width=6,height=2.8)
ggplot(s, aes(x = Reported_FDR, y = Real_FDR, color=method)) + 
  # geom_point(size=0.5)  +
  geom_line(size=0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = 'black')+
  # theme_tufte() +
    xlim(c(0,0.25))+
  # ylim(c(0,0.25))+
    # ylab('AUC (sensitivity-realFDR)') +
  scale_color_brewer(palette = 'Set1')  +
  facet_wrap(~comparison, nrow =2)+
  theme_classic()+
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust=1))
  
dev.off()



pdf(paste0(pdir, 'sensitivity_realFDR_strength4.pdf'),width=6,height=2.8)
ggplot(s, aes(x = Real_FDR, y = Sensitivity, color=method)) + 
  # geom_point(size=0.5)  +
  geom_line(size=0.5) + 
  # theme_tufte() +
    xlim(c(0,0.25))+
  # ylim(c(0,0.25))+
    # ylab('AUC (sensitivity-realFDR)') +
  scale_color_brewer(palette = 'Set1')  +
  facet_wrap(~comparison, nrow =2)+
  theme_classic()+
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust=1))
  
dev.off()


pdf(paste0(pdir, 'realFDR_reportedFDR_strength4_tendOnly.pdf'),width=4,height=2)
ggplot(s[s[,5]=='trendOnly', ], aes(x = Reported_FDR, y = Real_FDR,  color=method)) + 
  # geom_point(size=0.5)  +
  geom_line(size=0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = 'black')+
  # theme_tufte() +
    xlim(c(0,0.25))+
  # ylim(c(0,0.25))+
    # ylab('AUC (sensitivity-realFDR)') +
  scale_color_brewer(palette = 'Set1')  +
  # facet_wrap(~comparison, nrow =2)+
  theme_classic()+
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust=1))
dev.off()

pdf(paste0(pdir, 'sensitivity_realFDR_strength4_trendOnly.pdf'),width=4,height=2)
ggplot(s[s[,5]=='trendOnly', ], aes(x = Real_FDR, y = Sensitivity, color=method)) + 
  # geom_point(size=0.5)  +
  geom_line(size=0.5) + 
  # theme_tufte() +
    xlim(c(0,0.25))+
  # ylim(c(0,0.25))+
    # ylab('AUC (sensitivity-realFDR)') +
  scale_color_brewer(palette = 'Set1')  +
  # facet_wrap(~comparison, nrow =2)+
  theme_classic()+
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust=1))
dev.off()
