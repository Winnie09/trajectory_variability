source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/selgene/selgene.rds')
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/result/'
am <- list.files(ddir)


m = am[1]
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('testres', af)]


f = af[1]

df1 <- sapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = r[['res']]
  res = res[order(res[,1]), , drop=F]
  c(sub('.rds','',f), 'EM_SelectKnots', AreaUnderSensFdr(SensFdr(rownames(res), selgene, res)))
})

df1 = t(df1)
  
## tradeSeq
m = 'tradeSeq'
af = list.files(paste0(ddir, m, '/'))
af = af[grep('.rds', af)]
af = af[!grepl('sce', af)]
f = af[1]

df <- lapply(af, function(f){
  print(f)
  r = readRDS(paste0(ddir, m, '/', f))
  res = r[['startVsEndTest']]
  res = res[order(res[,'adj.P.Val']), , drop=F]
  a = c(sub('.rds','',f), 'tradeSeq_startVsEndTest', AreaUnderSensFdr(SensFdr(rownames(res), selgene, res)))
  res = r[['associationTest']]
  res = res[order(res[,'adj.P.Val']), , drop=F]
  b = c(sub('.rds','',f), 'tradeSeq_associationTest', AreaUnderSensFdr(SensFdr(rownames(res), selgene, res)))
  rbind(a, b)
})
df2 <- do.call(rbind, df)

res <- rbind(df1, df2)
colnames(res) <- c('Type', 'Method', 'Fdr.Diff', 'AUC')



#####
pd <- data.frame(res, SignalType = gsub('_.*','',sub('clusterType','',res[,1])), stringsAsFactors = F)
pd$SignalStreghth <- as.numeric(sapply(pd$Type, function(i) sub('.*_', '', i) ))

pd[,2] = as.factor(pd[,2])
pd[,5] = as.factor(pd[,5])
pd[,3] = as.numeric(pd[,3])
pd[,4] = as.numeric(pd[,4])
library(ggplot2)
library(gridExtra)
# pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/plot/ndPlainAll/compare_fdr_diff.pdf',width=8,height=4)
ggplot(pd, aes(x = SignalStreghth, y = Fdr.Diff, color=Method)) + geom_point(size=1)  + geom_line(size=0.1) + theme_classic() + ylab('fdr.diff(real~reported - 0.25*0.25/2)') +
  scale_color_brewer(palette = 'Set1') + facet_wrap(~SignalType, scales = 'free')
dev.off()
# pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/addsignal/plot/ndPlainAllcompare_auc.pdf',width=8,height=4)
ggplot(pd, aes(x = SignalStreghth, y = AUC, color=Method)) + geom_point(size=1)  + geom_line(size=0.1) + theme_classic() + ylab('Area Under Sensitivity Real_FDR curve') +
  scale_color_brewer(palette = 'Set1') + facet_wrap(~SignalType, scales = 'free')
dev.off()


#####
# af = c('clusterType1_1.rds', 'clusterType1_2.rds', 'clusterType1_3.rds', 'clusterType1_4.rds')
af = c('clusterType10_1.rds', 'clusterType10_2.rds', 'clusterType10_3.rds', 'clusterType10_4.rds')
par(mfrow=c(1,4))
for (f in af){
  r = readRDS(paste0(ddir, m, '/', f))
  res = r[['res']]
  res = res[order(res[,1]), , drop=F]
  a = SensFdr(rownames(res), selgene, res)
  plot(a[,2] ~ a[,3], xlab='Reported FDR', ylab = 'Real FDR', pch = 19, main=sub('.rds','',f), cex=.5)
}
  

##
Order <- rownames(res)
TruePositive = selgene
statistics = res
perf <- t(sapply(seq(1,length(Order)), function(i){
  num <- sum(Order[seq(1,i)] %in% selgene)
  c(num/length(TruePositive), (i - num)/i, statistics$adj.P.Val[i])
}))

a = perf
plot(a[,2] ~ a[,3], xlab='Reported FDR', ylab = 'Real FDR', pch = 19, main=sub('.rds','',f), cex=.5)

if (nrow(perf) > 1){
  for (i in seq(2, nrow(perf))){
    if (perf[i-1,2] > perf[i,2]) perf[i-1,2] <- perf[i,2]
  }
}
colnames(perf) <- c('Sensitivity','Real_FDR','Reported_FDR')
rbind(c(0,0,0),perf)

### plot the gene with largest fdr, but is a selgene
## plot this gene after saver
expr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/saver/',f))
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/pseudotime.rds')
pseudotime = pt[,2]
names(pseudotime) <- pt[,1]
design <- cbind(rep(1,8))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'
g = 'PDS5B:ENSG00000083642'
plotGene(testptObj = r, Gene = g, Mat = expr, Pseudotime = pseudotime, Cellanno = cellanno, Design = design,  Alpha=0.5, Size=0.1, PlotPoints = TRUE, FreeScale = FALSE, BySample = FALSE, type = 'Time')


## plot the original count, and added signal
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
mat <- readRDS('./testtime/data/data/null/hsc_mep_ery_saver.rds')
mat = 2^mat -1
clu <- readRDS('./testtime/data/data/null/geneCluster.rds')

clumat = mat[names(clu[clu==1]), ]
s <- names(sort(apply(clumat,1,sd)/rowMeans(clumat)))

j = 1
set.seed(12345)
addgene <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], length(selgene), replace = T)
thegene <- addgene[which(selgene == g)]
plotGene(testptObj = r, Gene = g, Mat = mat, Pseudotime = pseudotime, Cellanno = cellanno, Design = design,  Alpha=1, Size=0.5, PlotPoints = TRUE, FreeScale = FALSE, BySample = FALSE, type = 'Time')

thegene = "GSKIP:ENSG00000100744"
plotGene(testptObj = r, Gene = thegene, Mat = mat, Pseudotime = pseudotime, Cellanno = cellanno, Design = design,  Alpha=0.5, Size=0.1, PlotPoints = TRUE, FreeScale = FALSE, BySample = FALSE, type = 'Time')
    
