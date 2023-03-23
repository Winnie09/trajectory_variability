rm(list=ls())
library(here)
setwd(here())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
pdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
selgene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene.rds')
selgene1 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene1.rds') ## trendonly
selgene2 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene2.rds') ## meanonly
selgene3 <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene3.rds') ## trendMean

source('function/01_function.R')
slist <- list()
m = 'EM_pm'
sigstren = '1'
f = paste0(sigstren,'.rds')

print(f)
r = readRDS(paste0(ddir, m, '/', f))
res = data.frame(adj.P.Val = r$statistics[,'fdr.overall'], zscore = r$statistics[,'z.overall'], stringsAsFactors = F)
res = res[order(res[,1], res[,2]), ]
slist[['Lamian']] = rbind(cbind(SensFdr(selgene, res), method = 'Lamian', comparison = 'overall.test'),
                          cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'Lamian', comparison = 'trendOnly'),
                          cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'Lamian', comparison = 'meanOnly'),
                          cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'Lamian', comparison = 'trendMean'))

res = data.frame(adj.P.Val = r$statistics[,'fdr.trendDiff'], zscore = r$statistics[,'z.trendDiff'], stringsAsFactors = F)
res = res[order(res[,1], res[,2]), ]
res = res[complete.cases(res), ]
slist[['Lamian']] = rbind(slist[['Lamian']],cbind(SensFdr(c(selgene1,selgene3), res), method = 'Lamian', comparison = 'trend.test'))

res = data.frame(adj.P.Val = r$statistics[,'fdr.meanDiff'], zscore = r$statistics[,'z.meanDiff'], stringsAsFactors = F)
res = res[order(res[,1], res[,2]), ]
res = res[complete.cases(res), ]
slist[['Lamian']] = rbind(slist[['Lamian']],cbind(SensFdr(c(selgene2,selgene3), res), method = 'Lamian', comparison = 'mean.test'))

### Lamian.chisq
m = 'EM_chisq'
r = readRDS(paste0(ddir, m, '/', f))
res = data.frame(adj.P.Val = r$statistics[,'fdr.chisq.overall'], zscore = (r$ll3-r$ll1)/r$statistics[,'df.diff.overall'], stringsAsFactors = F)
res = res[order(res[,1], -res[,2]), ]
slist[['Lamian.chisq']] = rbind(cbind(SensFdr(selgene, res), method = 'Lamian.chisq', comparison = 'overall.test'),
                          cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'Lamian.chisq', comparison = 'trendOnly'),
                          cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'Lamian.chisq', comparison = 'meanOnly'),
                          cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'Lamian.chisq', comparison = 'trendMean'))

res = data.frame(adj.P.Val = r$statistics[,'fdr.chisq.trendDiff'], zscore = (r$ll3-r$ll2)/r$statistics[,'df.diff.trendDiff'], stringsAsFactors = F)
res = res[order(res[,1], -res[,2]), ]
res = res[complete.cases(res), ]
slist[['Lamian.chisq']] = rbind(slist[['Lamian.chisq']],cbind(SensFdr(c(selgene1,selgene3), res), method = 'Lamian.chisq', comparison = 'trend.test'))

res = data.frame(adj.P.Val = r$statistics[,'fdr.chisq.meanDiff'], zscore = (r$ll2 - r$ll1)/r$statistics[,'df.diff.meanDiff'], stringsAsFactors = F)
res = res[order(res[,1], -res[,2]), ]
res = res[complete.cases(res), ]
slist[['Lamian.chisq']] = rbind(slist[['Lamian.chisq']],cbind(SensFdr(c(selgene2,selgene3), res), method = 'Lamian.chisq', comparison = 'mean.test'))


###### phenopath100
fit = readRDS(paste0(ddir, 'phenopath100/fit_', sigstren, '.rds'))
zscore <- abs(fit$m_beta[1,]/sqrt(fit$s_beta[1,])) 
names(zscore) <- fit$feature_names
pval <- pnorm(zscore,lower.tail = F)
res = data.frame(score=zscore,pval=pval,fdr=p.adjust(pval,method='fdr'))
res <- res[order(-res[,1],res[,2]),]
slist[['phenopath']] = rbind(cbind(SensFdr(selgene, res), method = 'phenopath', comparison = 'overall.test'),
                             cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'phenopath', comparison = 'trendOnly'),
                             cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'phenopath', comparison = 'meanOnly'),
                             cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'phenopath', comparison = 'trendMean'))


## 'condiments'
res = readRDS(paste0(ddir, 'condiments/cond_genes_', sigstren, '.rds'))
res[is.na(res[,3]),3] <- 1
res$FDR <- p.adjust(res[,3],method='fdr')
res = res[order(res[, 3],-abs(res[, 1])),]
slist[['condiments']] = rbind(cbind(SensFdr(selgene, res), method = 'condiments', comparison = 'overall.test'),
                              cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'condiments', comparison = 'trendOnly'),
                              cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'condiments', comparison = 'meanOnly'),
                              cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'condiments', comparison = 'trendMean'))



## 'monocle2_trajtest'
res = readRDS(paste0(ddir, 'monocle2_trajtest/', sigstren, '.rds'))
res = res[order(res[, 1]),]
slist[['monocle2_trajTest']] = rbind(cbind(SensFdr(selgene, res), method = 'monocle2_trajTest', comparison = 'overall.test'),
                                     cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'monocle2_trajTest', comparison = 'trendOnly'),
                                     cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'monocle2_trajTest', comparison = 'meanOnly'),
                                     cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'monocle2_trajTest', comparison = 'trendMean'))

## 'monocle2_trajtest.corr'
res = readRDS(paste0(ddir, 'monocle2_trajtest.corr/', sigstren, '.rds'))
res = res[order(res[, 1]),]
slist[['monocle2_trajTest.corr']] = rbind(cbind(SensFdr(selgene, res), method = 'monocle2_trajTest.corr', comparison = 'overall.test'),
                                     cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'monocle2_trajTest.corr', comparison = 'trendOnly'),
                                     cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'monocle2_trajTest.corr', comparison = 'meanOnly'),
                                     cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'monocle2_trajTest.corr', comparison = 'trendMean'))


######## tradeSeq
m = 'tradeSeq'
res = readRDS(paste0(ddir, m, '/', f))[[1]]
res = res[order(res[, 3],-abs(res[, 1])),]
slist[['tradeSeq_diffEndTest']] = rbind(cbind(SensFdr(selgene, res), method = 'tradeSeq_diffEndTest', comparison = 'overall.test'),
                                        cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'tradeSeq_diffEndTest', comparison = 'trendOnly'),
                                        cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'tradeSeq_diffEndTest', comparison = 'meanOnly'),
                                        cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'tradeSeq_diffEndTest', comparison = 'trendMean'))

res = readRDS(paste0(ddir, m, '/', f))[[2]]
res[is.na(res[,3]),3] <- 1
res[is.na(res[,2]),2] <- 1
res[is.na(res[,1]),1] <- 0
res = res[order(res[, 3],-abs(res[, 1])),]
slist[['tradeSeq_patternTest']] = rbind(cbind(SensFdr(selgene, res), method = 'tradeSeq_patternTest', comparison = 'overall.test'),
                                        cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'tradeSeq_patternTest', comparison = 'trendOnly'),
                                        cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'tradeSeq_patternTest', comparison = 'meanOnly'),
                                        cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'tradeSeq_patternTest', comparison = 'trendMean'))

res = readRDS(paste0(ddir, m, '/', f))[[3]]
res[is.na(res[,3]),3] <- 1
res[is.na(res[,2]),2] <- 1
res[is.na(res[,1]),1] <- 0
res = res[order(res[, 3],-abs(res[, 1])),]
slist[['tradeSeq_earlyDETest']] = rbind(cbind(SensFdr(selgene, res), method = 'tradeSeq_earlyDETest', comparison = 'overall.test'),
                                        cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'tradeSeq_earlyDETest', comparison = 'trendOnly'),
                                        cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'tradeSeq_earlyDETest', comparison = 'meanOnly'),
                                        cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'tradeSeq_earlyDETest', comparison = 'trendMean'))

## limma
m = 'meandiff'
res = readRDS(paste0(ddir, m, '/', f))
res = res[order(res[,5], -abs(res[,1])), ]
slist[['limma']] = rbind(cbind(SensFdr(selgene, res), method = 'limma', comparison = 'overall.test'),
                         cbind(SensFdr(selgene1, res[!rownames(res) %in% c(selgene2, selgene3), ]), method = 'limma', comparison = 'trendOnly'),
                         cbind(SensFdr(selgene2, res[!rownames(res) %in% c(selgene1, selgene3), ]), method = 'limma', comparison = 'meanOnly'),
                         cbind(SensFdr(selgene3, res[!rownames(res) %in% c(selgene1, selgene2), ]), method = 'limma', comparison = 'trendMean'))


## concatenate
saveRDS(slist, paste0(rdir,'sensfdr_list_', sigstren, '.rds'))
s <- do.call(rbind, slist)
str(s)
s <- as.data.frame(s)
s[,1] <- as.numeric(as.character(s[,1]))
s[,3] <- as.numeric(as.character(s[,3]))
s[,2] <- as.numeric(as.character(s[,2]))
s[,5] <- factor(as.character(s[,5]), levels = c("overall.test", "trend.test","mean.test","trendOnly","trendMean", "meanOnly" ))
saveRDS(s, paste0(rdir,'sensfdr_', sigstren, '.rds'))


##########
## plot ##
##########
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
a = read.csv(paste0(pdir, 'color_code.csv'), row.names = 1, stringsAsFactors = F)
colv <- a[,2]
names(colv) <- a[,1]
 
pdf(paste0(pdir, 'realFDR_reportedFDR_strength_', sigstren, '.pdf'),width=6,height=3.5)
ggplot(s, aes(x = Reported_FDR, y = Real_FDR, color=method)) + 
  geom_line(size=0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = 'black')+
  xlim(c(0,0.25))+
  facet_wrap(~comparison, nrow =2)+
  scale_color_manual(values = colv)  + theme_compact() +
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust=0))
dev.off()

pdf(paste0(pdir, 'sensitivity_realFDR_strength', sigstren, '.pdf'),width=6,height=3.5)
ggplot(s, aes(x = Real_FDR, y = Sensitivity, color=method)) + 
  geom_line(size=0.5) + 
  xlim(c(0,0.25))+
  facet_wrap(~comparison, nrow =2)+
  scale_color_manual(values = colv)  + theme_compact() +
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust=0))

dev.off()

pdf(paste0(pdir, 'realFDR_reportedFDR_strength_', sigstren, '_tendOnly.pdf'),width=3.9,height=2.2)
ggplot(s[s[,5]=='trendOnly', ], aes(x = Reported_FDR, y = Real_FDR,  color=method)) + 
  geom_line(size=0.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = 'black')+
  xlim(c(0,0.25))+
  scale_color_manual(values = colv)  + theme_compact() +
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust = 0))
dev.off()

pdf(paste0(pdir, 'sensitivity_realFDR_strength_', sigstren, '_trendOnly.pdf'),width=3.9,height=2.2)
ggplot(s[s[,5]=='trendOnly', ], aes(x = Real_FDR, y = Sensitivity, color=method)) + 
  geom_line(size=0.5) + 
  xlim(c(0,0.25))+
  scale_color_manual(values = colv)  + theme_compact() +
  theme(legend.position = 'right', axis.text.x = element_text(angle=45, hjust = 0))
dev.off()

