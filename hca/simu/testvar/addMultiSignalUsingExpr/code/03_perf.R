rm(list=ls())
library(here)
setwd(here())
#setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- 'hca/simu/testvar/addMultiSignalUsingExpr/result/perf/'
selgene <- readRDS('hca/simu/testvar/addMultiSignalUsingExpr/data/selgene/selgene.rds')
source('function/01_function.R')

perflist <- list()
m = 'EM_pm'
af = list.files(paste0(ddir, m, '/'))
af <- af[grepl('.rds', af)]
df1.nc <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  
  res = data.frame(adj.P.Val = r$statistics[,'fdr.overall'], zscore = r$statistics[,'z.overall'], stringsAsFactors = F)
  res = res[order(res[,1], res[,2]), ]
  c(sub('.rds','',f), 'Lamian', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['Lamian']] = t(df1.nc)

sav <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/data/saver/0.5.rds')
m = 'phenopath100'
af = list.files(paste0(ddir, m, '/'),pattern = 'fit')
# use z-score to order, same as in its function significant_interactions().
df3 <- sapply(af, function(f) {
  fit = readRDS(paste0(ddir, m, '/', f))
  sig = readRDS(paste0(ddir, m, '/', sub('fit','sig',f)))[,1]
  names(sig) = rownames(sav)
  zscore <- abs(fit$m_beta[1,]/sqrt(fit$s_beta[1,])) 
  names(zscore) <- fit$feature_names
  pval <- pnorm(zscore,lower.tail = F)
  res = data.frame(score=zscore,pval=pval,fdr=p.adjust(pval,method='fdr'))
  res <- res[order(-res[,1],res[,2]),]
  a = setdiff(names(sig)[sig], rownames(res)[1:1499])
  #  > str(a)
  # chr(0) 
  c(sub('fit_','',sub('.rds', '', f)), 'phenopath', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['phenopath']] = t(df3)

m = 'condiments'
af = list.files(paste0(ddir, m, '/'),pattern = 'cond_genes')
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))
  res[is.na(res[,3]),3] <- 1
  res$FDR <- p.adjust(res[,3],method='fdr')
  res = res[order(res[, 3],-abs(res[, 1])),]
  c(sub('cond_genes_','',sub('.rds', '', f)), 'condiments', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['condiments']] = t(df3)

m = 'monocle2_trajtest'
af = list.files(paste0(ddir, m, '/'))
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[, 1]),]
  c(sub('.rds', '', f), 'monocle2_trajTest', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['monocle2_trajTest']] = t(df3)


m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
af = af[!grepl('sce', af)]
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))[[1]]
  res = res[order(res[, 3],-abs(res[, 1])),]
  c(sub('.rds', '', f), 'tradeSeq_diffEndTest', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tradeSeq_diffEndTest']] = t(df3)

m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
af = af[!grepl('sce', af)]
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))[[2]]
  res[is.na(res[,3]),3] <- 1
  res[is.na(res[,2]),2] <- 1
  res[is.na(res[,1]),1] <- 0
  res = res[order(res[, 3],-abs(res[, 1])),]
  c(sub('.rds', '', f), 'tradeSeq_patternTest', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tradeSeq_patternTest']] = t(df3)

m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
af = af[!grepl('sce', af)]
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))[[3]]
  res[is.na(res[,3]),3] <- 1
  res[is.na(res[,2]),2] <- 1
  res[is.na(res[,1]),1] <- 0
  res = res[order(res[, 3],-abs(res[, 1])),]
  c(sub('.rds', '', f), 'tradeSeq_earlyDETest', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['tradeSeq_earlyDETest']] = t(df3)


m = 'meandiff'
af = list.files(paste0(ddir, m, '/'))
af <- af[grepl('.rds', af)]
df2 <- sapply(af, function(f){
  print(f)
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[,5], -abs(res[,1])), ]
  c(sub('.rds','',f), 'limma', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['limma']] = t(df2)

m = 'EM_chisq'
af = list.files(paste0(ddir, m, '/'))
af <- af[grepl('.rds', af)]
df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = data.frame(fdr = r$statistics[,'fdr.chisq.overall'], stat = (r$ll3-r$ll1)/r$statistics[,'df.diff.overall'], stringsAsFactors = F)
  res = res[order(res[,1], -res[,2]), ]
  c(sub('.rds','',f), 'Lamian.chisq', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['Lamian_chisq']] = t(df1)

m = 'monocle2_trajtest.corr'
af = list.files(paste0(ddir, m, '/'))
df3 <- sapply(af, function(f) {
  res = readRDS(paste0(ddir, m, '/', f))
  res = res[order(res[, 1]),]
  c(sub('.rds', '', f), 'monocle2_trajTest.corr', AreaUnderSensFdr(SensFdr(selgene, res)))
})
perflist[['monocle2_trajTest.corr']] = t(df3)

## concatenate
saveRDS(perflist, paste0(rdir,'perflist.rds'))
perf <- do.call(rbind, perflist)
colnames(perf) <- c('SignalStrength', 'Method', 'Fdr.Diff', 'AUC')
perf <- as.data.frame(perf)
perf[,1] <- as.numeric(as.character(perf[,1]))
perf[,3] <- as.numeric(as.character(perf[,3]))
perf[,4] <- as.numeric(as.character(perf[,4]))
saveRDS(perf, paste0(rdir,'perf.rds'))
rm(list=ls())


