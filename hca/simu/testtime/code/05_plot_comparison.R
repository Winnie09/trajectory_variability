source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/selgene/selgene.rds')
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/result/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/result/perf/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/plot/'
res = readRDS(paste0(rdir, 'perf.rds'))

### plot
pd <- data.frame(res, SignalType = gsub('_.*','',sub('clusterType','',res[,1])), stringsAsFactors = F)
pd$SignalStreghth <- as.numeric(sapply(pd$Type, function(i) sub('.*_', '', i) ))
pd[,2] = as.factor(pd[,2])
pd[,5] = as.factor(pd[,5])
pd[,3] = as.numeric(pd[,3])
pd[,4] = as.numeric(pd[,4])
library(ggplot2)
library(gridExtra)
pdf(paste0(pdir, 'compare_fdr_diff.pdf'),width=8,height=4)
ggplot(pd, aes(x = SignalStreghth, y = Fdr.Diff, color=Method)) + 
  geom_point(size=1)  + 
  geom_line(size=0.1) + 
  theme_classic() + 
  ylab('fdr.diff(real~reported - 0.25*0.25/2)') +
  facet_wrap(~SignalType, scales = 'free')
dev.off()
pdf(paste0(pdir, 'compare_auc.pdf'), width=8,height=4)
ggplot(pd, aes(x = SignalStreghth, y = AUC, color=Method)) + 
  geom_point(size=1)  + 
  geom_line(size=0.1) + 
  theme_classic() + ylab('Area Under Sensitivity Real_FDR curve') +
  scale_color_brewer(palette = 'Set1') + 
  facet_wrap(~SignalType, scales = 'free')
dev.off()


### plot curves
af = c('clusterType1_1.rds', 'clusterType1_2.rds', 'clusterType1_3.rds', 'clusterType1_4.rds')
# af = c('clusterType10_1.rds', 'clusterType10_2.rds', 'clusterType10_3.rds', 'clusterType10_4.rds')
par(mfrow=c(2,4))
m = 'EM_SelectKnots'
for (f in af){
  r = readRDS(paste0(ddir, m, '/', f))
  res = r[['res']]
  fdrchar <- intersect(colnames(res), c('adj.P.Val','adj.pvalue','fdr','FDR','Fdr','adj.p', 'adj.P','adj.Pval'))
  fdrcol <- which(colnames(res) == fdrchar)
  res = res[order(res[,fdrcol]), , drop=F]
  a = SensFdr(selgene, res)
  plot(a[,2] ~ a[,3], xlab='Reported FDR', ylab = 'Real FDR', pch = 19, main=sub('.rds','',f), cex=.5)
}
  
for (f in af){
  r = readRDS(paste0(ddir, m, '/', f))
  res = r[['res']]
  fdrchar <- intersect(colnames(res), c('adj.P.Val','adj.pvalue','fdr','FDR','Fdr','adj.p', 'adj.P','adj.Pval'))
  fdrcol <- which(colnames(res) == fdrchar)
  res = res[order(res[,fdrcol]), , drop=F]
  a = SensFdr(selgene, res)
  plot(a[,1] ~ a[,2], xlab='Real FDR', ylab = 'Sensitivity', pch = 19, main=sub('.rds','',f), cex=.5)
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


