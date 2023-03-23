library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(here)
library(parallel)
setwd(here())
compareid <- (1:50)*100
scompareid <- (1:2)*100
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/res/'
# af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/res',full.names = T)
# d1 <- lapply(af,readRDS)
# expid <- expand.grid(1:length(d1),1:length(d1))
# expid <- expid[expid[,1] < expid[,2],]
# ao1 <- sapply(1:nrow(expid),function(i) {
#   s1 <- data.frame(d1[[unlist(expid[i,1])]][[1]])
#   s2 <- data.frame(d1[[unlist(expid[i,2])]][[1]])
#   s1 <- s1[order(s1$pval.overall,-s1$z.overall),]
#   s2 <- s2[order(s2$pval.overall,-s2$z.overall),]
#   mean(sapply(scompareid,function(i) {
#     length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
#   }))
# })
# i <- which.max(ao1)
# s1 <- data.frame(d1[[unlist(expid[i,1])]][[1]])
# s2 <- data.frame(d1[[unlist(expid[i,2])]][[1]])
# s1 <- s1[order(s1$pval.overall,-s1$z.overall),]
# s2 <- s2[order(s2$pval.overall,-s2$z.overall),]
# saveRDS(list(s1, s2), paste0(rdir, 'Lamian.rds'))
res <- readRDS(paste0(rdir, 'Lamian.rds'))
Lamian.pm.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

############### Lamian.chisq
# af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/compres/Lamian.chisq/',full.names = T)
# d1 <- mclapply(af,function(f) {
#     k <- readRDS(f)
#     data.frame(k[[1]], stat = (k$ll3-k$ll1)/k$statistics[,'df.diff.overall'])
# },mc.cores=10)
# expid <- expand.grid(1:length(d1),1:length(d1))
# expid <- expid[expid[,1] < expid[,2],]
# ao1 <- sapply(1:nrow(expid),function(i) {
#   s1 = d1[[unlist(expid[i,1])]]
#   s2 = d1[[unlist(expid[i,2])]]
#   s1 <- s1[order(s1$pval.chisq.overall,-s1$stat),] #######
#   s2 <- s2[order(s2$pval.chisq.overall,-s2$stat),]
#   sapply(compareid,function(i) {
#     length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
#   })
# })
# i <- apply(ao1,2,function(i) sum(abs(i-Lamian.pm.o)))
# i <- which(i==sort(i)[1])
# s1 = d1[[unlist(expid[i,1])]]
# s2 = d1[[unlist(expid[i,2])]]
# s1 <- s1[order(s1$pval.chisq.overall,-s1$stat),] #######
# s2 <- s2[order(s2$pval.chisq.overall,-s2$stat),]
# saveRDS(list(s1, s2), paste0(rdir, 'Lamian.chisq.rds'))
res <- readRDS(paste0(rdir, 'Lamian.chisq.rds'))
Lamian.chisq.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

############ Lamian_TDE.pm
# af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/compres/Lamian_TDE.pm/',full.names = T, pattern = '.rds$')
# d1 <- mclapply(af,function(i) {
#   s1 <- data.frame(readRDS(i)[[1]])
#   s1[order(s1$pval.overall,-s1$z.overall),]
# },mc.cores=3)
# expid <- expand.grid(1:length(d1),1:length(d1))
# expid <- expid[expid[,1] < expid[,2],]
# # ao1 <- sapply(1:nrow(expid),function(i) {
# #   s1 <- d1[[unlist(expid[i,1])]]
# #   s2 <- d1[[unlist(expid[i,2])]]
# #   mean(sapply(compareid,function(i) {
# #     length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
# #   }))
# # })
# 
# ao1 <- sapply(1:nrow(expid),function(i) {
#   s1 <- d1[[unlist(expid[i,1])]]
#   s2 <- d1[[unlist(expid[i,2])]]
#   sapply(compareid,function(i) {
#     length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
#   })
# })
# 
# i <- which.max(colSums(ao1[1:4,]))
# s1 <- d1[[unlist(expid[i,1])]]
# s2 <- d1[[unlist(expid[i,2])]]
# saveRDS(list(s1, s2), paste0(rdir, 'Lamian_TDE.pm.rds'))
res <- readRDS(paste0(rdir, 'Lamian_TDE.pm.rds'))
Lamian_tde.pm.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

############### Lamian_TDE.chisq
# af <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/compres/Lamian_TDE.chisq',full.names = T, pattern = '.rds$')
# d1 <- lapply(af,readRDS)
# expid <- expand.grid(1:length(d1),1:length(d1))
# expid <- expid[expid[,1] < expid[,2],]
# ao1 <- sapply(1:nrow(expid),function(i) {
#   s1 = d1[[unlist(expid[i,1])]][[1]]
#   s2 = d1[[unlist(expid[i,2])]][[1]]
#   s1 <- s1[order(s1$pval.chisq.overall,-s1$llr),] #######
#   s2 <- s2[order(s2$pval.chisq.overall,-s2$llr),]
#   mean(sapply(scompareid,function(i) {
#     length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
#   }))
# })
# i <- which.max(ao1)
# s1 = d1[[unlist(expid[i,1])]][[1]]
# s2 = d1[[unlist(expid[i,2])]][[1]]
# s1 <- s1[order(s1$pval.chisq.overall,-s1$llr),] #######
# s2 <- s2[order(s2$pval.chisq.overall,-s2$llr),]
# saveRDS(list(s1, s2), paste0(rdir, 'Lamian.TDE.chisq.rds'))
res <- readRDS(paste0(rdir, 'Lamian.TDE.chisq.rds'))
Lamian_tde.chisq.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})


############## tradeseq
# af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/tradeseq',full.names = T)
# d2 <- lapply(af,readRDS)
# expid <- expand.grid(1:length(d2),1:length(d2))
# expid <- expid[expid[,1] < expid[,2],]
# ao2 <- sapply(1:3,function(metid) {
#   sapply(1:nrow(expid),function(i) {
#     s1 <- data.frame(d2[[unlist(expid[i,1])]][[metid]])
#     s2 <- data.frame(d2[[unlist(expid[i,2])]][[metid]])
#     mean(sapply(scompareid,function(i) {
#       length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
#     }))
#   }) 
# })
# i <- which.min(apply(ao2,1,max))
# for (metid in 1:3) {
#   s1 <- data.frame(d2[[unlist(expid[i,1])]][[metid]])
#   s2 <- data.frame(d2[[unlist(expid[i,2])]][[metid]])
#   saveRDS(list(s1, s2), paste0(rdir, 'tradeSeq_', names(d2[[unlist(expid[i,1])]])[metid],'.rds'))
# }
tradeseq.o <- sapply(list.files(rdir,pattern = 'tradeSeq'),function(metid) {
  res <- readRDS(paste0(rdir, metid)) ########
  o <- sapply(compareid,function(i) {
    length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
  })
})

##############
# af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/limma',full.names = T)[1:30]
# ao3 <- sapply(af,function(i) {
#   res <- readRDS(i)
#   s1 <- res[[1]]
#   s2 <- res[[2]]
#   mean(sapply(scompareid,function(i) {
#     length(intersect(rownames(s1)[1:i],rownames(s2)[1:i]))/i
#   }))
# }) 
# i <- af[which.min(ao3)]
# res <- readRDS(i)
# s1 <- res[[1]]
# s2 <- res[[2]]
# saveRDS(list(s1, s2), paste0(rdir, 'limma.rds'))
res <- readRDS(paste0(rdir, 'limma.rds'))
limma.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

###################
# af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/condiments',full.names = T, pattern = 'cond_gene')
# af <- af[1:min(30, length(af))]
# af = unique(sapply(af, function(i) strsplit(i, '_')[[1]][1]))
# ao4 <- sapply(af,function(i) {
#   s1 <- readRDS(paste0(i, '_1cond_gene_res.rds'))
#   s1[is.na(s1[,3]),3] <- 1
#   s1$FDR <- p.adjust(s1[,3],method='fdr')
#   s1 = s1[order(s1[, 3],-abs(s1[, 1])),]
#   s2 <- readRDS(paste0(i, '_2cond_gene_res.rds'))
#   s2[is.na(s2[,3]),3] <- 1
#   s2$FDR <- p.adjust(s2[,3],method='fdr')
#   s2 = s2[order(s2[, 3],-abs(s2[, 1])),]
#   mean(sapply(scompareid,function(j) {
#     length(intersect(rownames(s1)[1:j],rownames(s2)[1:j]))/j
#   }))
# }) 
# i <- af[which.min(ao4)]
# s1 <- readRDS(paste0(i, '_1cond_gene_res.rds'))
# s1[is.na(s1[,3]),3] <- 1
# s1$FDR <- p.adjust(s1[,3],method='fdr')
# s1 = s1[order(s1[, 3],-abs(s1[, 1])),]
# s2 <- readRDS(paste0(i, '_2cond_gene_res.rds'))
# s2[is.na(s2[,3]),3] <- 1
# s2$FDR <- p.adjust(s2[,3],method='fdr')
# s2 = s2[order(s2[, 3],-abs(s2[, 1])),]
# saveRDS(list(s1, s2), paste0(rdir, 'condiments.rds'))
res <- readRDS(paste0(rdir, 'condiments.rds'))
condiments.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

###  monocle2 trajtest
# af <- list.files('/home-4/zji4@jhu.edu/scratch/diffpt/su/repro/compres/monocle2_trajtest',full.names = T)
# af <- af[1:min(30, length(af))]
# af = unique(sapply(af, function(i) paste0(strsplit(i, '_')[[1]][1], '_', strsplit(i, '_')[[1]][2])))
# ao5 <- sapply(af,function(i) {
#   s1 <- readRDS(paste0(i, '_1.rds'))
#   s1 = s1[order(s1[,3], -s1[, 1]),]
#   s2 <- readRDS(paste0(i, '_2.rds'))
#   s2 = s2[order(s2[,3], -s2[, 1]),]
#   mean(sapply(compareid,function(j) {
#     length(intersect(rownames(s1)[1:j],rownames(s2)[1:j]))/j
#   }))
# })
# i <- af[which.min(ao5)]
# s1 <- readRDS(paste0(i, '_1.rds'))
# s1 = s1[order(s1[,3], -s1[, 1]),]
# s2 <- readRDS(paste0(i, '_2.rds'))
# s2 = s2[order(s2[,3], -s2[, 1]),]
# saveRDS(list(s1, s2), paste0(rdir, 'monocle2_trajTest.rds'))
res <- readRDS(paste0(rdir, 'monocle2_trajTest.rds'))
monocle2_trajTest.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

### monocle2 trajtest correction
# ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/compres/monocle2_trajtest.corr/'
# af <- list.files(ddir, full.names = T)
# af <- af[1:min(30, length(af))]
# af = unique(sapply(af, function(i) sub('_.*', '',sub('.*/','', i))))
# ao5 <- sapply(af,function(i) {
#   s1 <- readRDS(paste0(ddir, i, '_1.rds'))
#   s1 = s1[order(s1[,3], -s1[, 1]),]
#   s2 <- readRDS(paste0(ddir, i, '_2.rds'))
#   s2 = s2[order(s2[,3], -s2[, 1]),]
#   mean(sapply(scompareid,function(j) {
#     length(intersect(rownames(s1)[1:j],rownames(s2)[1:j]))/j
#   }))
# }) 
# i <- af[which.min(ao5)]
# s1 <- readRDS(paste0(ddir, i, '_1.rds'))
# s1 = s1[order(s1[,3], -s1[, 1]),]
# s2 <- readRDS(paste0(ddir, i, '_2.rds'))
# s2 = s2[order(s2[,3], -s2[, 1]),]
# saveRDS(list(s1, s2), paste0(rdir, 'monocle2_trajTest.corr.rds'))
res <- readRDS(paste0(rdir, 'monocle2_trajTest.corr.rds'))
monocle2_trajTest.corr.o <- sapply(compareid,function(i) {
  length(intersect(rownames(res[[1]])[1:i],rownames(res[[2]])[1:i]))/i
})

saveRDS(list(Lamian.pm=Lamian.pm.o,Lamian.chisq = Lamian.chisq.o, Lamian_tde.pm = Lamian_tde.pm.o, Lamian_tde.chisq = Lamian_tde.chisq.o, tradeseq=tradeseq.o, limma=limma.o, condiments=condiments.o, monocle2_trajTest=monocle2_trajTest.o, monocle2_trajTest.corr=monocle2_trajTest.corr.o), paste0(rdir, 'overlap.rds'))
####### ======== phenopath10, try the following code when ready
# fit <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/compres/phenopath10/fit_res.rds'))
# zscore <- abs(fit$m_beta[1,]/sqrt(fit$s_beta[1,]))
# names(zscore) <- fit$feature_names
# pval <- pnorm(zscore,lower.tail = F)
# res = data.frame(score=zscore,pval=pval,fdr=p.adjust(pval,method='fdr'))
# res <- res[order(-res[,1],res[,2]),]
# phenopath_genes = rownames(res)
# phenopath_sig = rownames(res)[res[,'fdr'] < 0.05]
##################3
df <- cbind(Lamian.pm.o, Lamian.chisq.o, Lamian_tde.pm.o,Lamian_tde.chisq.o, tradeseq.o, limma.o, condiments.o,monocle2_trajTest.o, monocle2_trajTest.corr.o)
colnames(df) <- c('Lamian.pm','Lamian.chisq','Lamian_TDE.pm','Lamian_TDE.chisq', sub('.rds','',colnames(tradeseq.o)),'Limma', 'condiments', 'monocle2_trajTest', 'monocle2_trajTest.corr')
df <- melt(df)
colnames(df) <- c('numgene','method','overlap')
df[,1] <- compareid
saveRDS(df,'/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/df.rds')

df = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/df.rds')
mycol = brewer.pal(12, 'Paired')[1:length(unique(df[,2]))]
names(mycol) = unique(as.character(df[,2]))

# pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/overlapline.pdf',width=4,height=2.1)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/overlapline.pdf',width=5,height=3.3)
ggplot(df,aes(x=numgene,y=overlap,col=method)) + 
  geom_point(size = 0.2) + 
  geom_line(size = 0.1) + theme_classic() +
  xlab('number of top genes') + ylab('overlap proportion') +
  scale_color_manual(values = mycol) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pd <- tapply(df[,3],df[,2],mean)
pd <- data.frame(method=names(pd),overlap=pd)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/repro/plot/overlap.pdf',width=4,height=4)
ggplot(pd,aes(method,overlap)) + geom_bar(stat='identity') + theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  ylab('average overlap')
dev.off()




