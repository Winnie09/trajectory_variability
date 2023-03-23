clusterType = 9
pctGene = 4
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('./function/01_function.R')
ddir <- './hca/data/simu/testtime/'
rdir <- './hca/simu/testtime/result/'
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))
testres = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/result/addsignal/EM_centered/clusterType', clusterType, '_', pctGene, '.rds'))
oldres <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/old_without_demean/result/addsignal/EM_SelectKnots/clusterType', clusterType, '_', pctGene, '_testres.rds'))
pseudotime <- readRDS(paste0(ddir, 'null/pseudotime.rds'))
expr <- readRDS(paste0(ddir, 'saver/clusterType', clusterType, '_', pctGene, '.rds'))


## one group along pseudotime
# ------------------------------------------------
# see if  demean the data will change the results
identical(names(testres$fdr), names(oldres$fdr))

# compare fdr
plot(oldres$fdr ~ testres$fdr, pch = 20)
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/selgene/selgene.rds')

pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testtime/plot/compare_demean_and_not_demean_clusterType', clusterType, '_', pctGene, '.pdf'), width = 5, height = 3.5)
plot(oldres$fdr ~ testres$fdr, pch = 20, col = ifelse(names(oldres$fdr) %in% selgene, 'red', 'black'), xlab = 'fdr (centered data)', ylab = 'fdr (un-centered data)')
dev.off()

# plot example genes
g <- names(which(testres$fdr < 0.05 & oldres$fdr > 0.5))
plotGene(testres, g[g %in% selgene])
plotGene(testres, g[!g %in% selgene])

g <- names(which(testres$fdr > 0.2 & oldres$fdr < 0.05))
g <- names(which(testres$fdr > 0.63 & oldres$fdr < 0.6))
g <- names(which(testres$fdr < 1e-15 & oldres$fdr > 0.5))
g <- names(which(testres$fdr < 1e-15 & oldres$fdr <1e-15))
expr.ori <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testtime/null/hsc_mep_ery_saver.rds')
plotGene(testres, g[g %in% selgene])

length(names(which(testres$fdr < 0.2 & oldres$fdr < 0.2)))
length(names(which(testres$fdr < 0.05 & oldres$fdr < 0.05)))

AreaUnderSensFdr(SensFdr(selgene, data.frame(fdr = testres$fdr, foldchange = testres$foldchange)))

# compare log-likelihood
expr.ori = testres$expr.ori
cellanno=testres$cellanno
design=testres$design
pseudotime=testres$pseudotime
rowm <- sapply(unique(cellanno[,2]), function(s){
  tmp = expr.ori[selgene, cellanno[,2] == s]
  rowMeans(tmp)
})

ag = rownames(rowm)[testres$fdr[rownames(rowm)] > 0.05]
min = sapply(ag, function(g) min(rowm[g,]))
min.all = sapply(1:nrow(rowm), function(g) min(rowm[g,]))


ag = rownames(expr.ori)[testres$fdr > 0.05 & oldres$fdr < 0.05]
ag = ag[ag %in% selgene]
summary(testres$orill[ag] - oldres$orill[ag])

g =  "DHX37:ENSG00000150990"
g = names(testres$fdr)[grep('^RPL39:', names(testres$fdr))]
par(mfrow=c(2,2))
plot(density(oldres$perll[g,]), main = 'EM_NOT_Centered', xlim=c(15400, 16500))
abline(v = oldres$orill[g], col='red')

plot(density(testres$perll[g,]), main = 'EM_Centered')
abline(v = testres$orill[g], col='red')

# -------------------------------------------------------
# plot orignal signals, added signals, and final signals
# -------------------------------------------------------
geneProp <- 0.2
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './hca/data/simu/testtime/'
### load saver, and count matrix
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
saverlog <- readRDS('./hca/data/simu/testtime/null/hsc_mep_ery_saver.rds')
cnt <- readRDS('./hca/data/simu/testtime/null/hsc_mep_ery_count.rds')

### prepare count, imputed, selected genes
savercnt <- 2^saverlog - 1
rownames(cnt) <- sapply(rownames(cnt), function(i) sub('_','-',i))
rownames(savercnt) <- sapply(rownames(savercnt), function(i) sub('_','-',i))
cnt <- cnt[, colnames(savercnt)]
allp = sub(':.*','', colnames(savercnt))
sample <- sub(':.*','',colnames(savercnt))
names(sample) <- colnames(savercnt)
set.seed(12345)
selgene <- sample(row.names(savercnt), round(geneProp * nrow(savercnt)))
othgene <- setdiff(rownames(savercnt), selgene)

### permute the pseudotime for each sample
pmlist <- lapply(unique(allp), function(p){
  tmp <- savercnt[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmsavercnt <- do.call(cbind, pmlist)

pmlist <- lapply(unique(allp), function(p){ ###### -->
  tmp <- cnt[,allp==p]
  set.seed(12345)
  colnames(tmp) <- sample(colnames(tmp))
  tmp
})
pmcnt <- do.call(cbind, pmlist)
identical(colnames(pmsavercnt), colnames(pmcnt)) ######### <------
identical(rownames(pmsavercnt), rownames(pmcnt)) ######### <------

# ========  
### order cells by pseudotime, and get the optimal clusters using saver imputed 80% genes
pt <- readRDS('./hca/data/simu/testtime/null/pseudotime.rds')
savercnt <- savercnt[, pt[,1]]
tmp <- savercnt[othgene, pt[,1]]

library(matrixStats)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- rowSds(data, center = cm)
  (data - cm) / csd
}

set.seed(12345)
tmp <- scalematrix(tmp)
library(parallel)
clu = readRDS(paste0(rdir, 'null/geneCluster.rds'))

### add signal to permuted expression (both cnt and savercnt) 
i = clusterType
j = pctGene
clumat = savercnt[names(clu[clu==i]), ]
# s <- names(sort(apply(clumat,1,sd)/rowMeans(clumat))) ####
sds <- sapply(unique(allp), function(p){
  tmp <- clumat[, allp == p]
  sd <- apply(tmp, 1, sd)   
  rank(sd)
})
s <- names(sort(rowMeans(sds), decreasing = TRUE))

# j = 1:4
set.seed(12345)
addgene <- sample(s[(1+0.25*(j-1)*length(s)) : (0.25*j*length(s))], length(selgene), replace = T)
addsavercnt <- savercnt[addgene, pt[,1]]
ressavercnt <-  pmsavercnt[selgene, pt[,1]] + addsavercnt
mat <- rbind(ressavercnt[selgene, pt[,1]], pmsavercnt[othgene, pt[,1]])
# saveRDS(mat, paste0(rdir, 'saver/clusterType', i, '_', j,'.rds'))
# saveRDS(list(selgene = selgene, addgene = rownames(addsavercnt)), paste0(rdir, 'saver/clusterType', i, '_', j,'_selgene_addgene.rds'))

addsavercnt <- cnt[addgene, pt[,1]]
ressavercnt <-  pmcnt[selgene, pt[,1]] + addsavercnt
mat <- rbind(ressavercnt[selgene, pt[,1]], pmcnt[othgene, pt[,1]])
# saveRDS(mat, paste0(rdir, 'count/clusterType', i, '_', j,'.rds'))

# plot example genes
g <- "RPL39:ENSG00000198918" 
g <- rownames(cnt)[grepl('^DHX37:', rownames(cnt))]
par(mfrow=c(3,3))
plot(addsavercnt[which(selgene == g), ] ~ seq(1, ncol(addsavercnt)), pch = 20, cex=.1)
plot(pmsavercnt[g,] ~ seq(1, ncol(pmsavercnt)), pch = 20, cex=.1)
plot(ressavercnt[g, ] ~ seq(1, ncol(ressavercnt)), pch = 20, cex=.1)

plot(addsavercnt[which(selgene == g), pt[,1]] ~ seq(1, ncol(addsavercnt)), pch = 20, cex=.1)
plot(pmsavercnt[g, pt[,1]] ~ seq(1, ncol(pmsavercnt)), pch = 20, cex=.1)
plot(ressavercnt[g, pt[,1]] ~ seq(1, ncol(ressavercnt)), pch = 20, cex=.1, col = allp)


library(ggplot2)
library(reshape2)
pd1 = data.frame(cnt = pmsavercnt[g, pt[,1]], 
                pseudotime = pt[,2], type = 'pmsavercnt', sample = sub(':.*', '', pt[,1]))
pd2 = data.frame(cnt = addsavercnt[which(selgene == g), pt[,1]],
                pseudotime = pt[,2], type = 'addsavercnt', sample = sub(':.*', '', pt[,1]))
pd3 = data.frame(cnt = ressavercnt[g, pt[,1]],
                pseudotime = pt[,2], type = 'ressavercnt', sample = sub(':.*', '', pt[,1]))
ggplot(data = rbind(pd1, pd2, pd3)) + geom_point(aes(y = cnt, x = pseudotime, color = sample), size = 0.1) +
  theme_classic() +
  facet_wrap(~type) +
  xlab('pseudotime') + ylab('savercount') +
  theme(axis.text.x = element_blank())

pd = rbind(pd1, pd2, pd3)
pd = pd[order(pd$sample, pd$pseudotime), ]
pd$order = rep(1:nrow(pd1), each = 3)
ggplot(data = pd) + geom_point(aes(y = cnt, x = order, color = sample, group = sample), size = 0.1) +
  theme_classic() +
  facet_wrap(~type) +
  xlab('pseudotime(stratify sample)') + ylab('savercount') +
  theme(axis.text.x = element_blank())


