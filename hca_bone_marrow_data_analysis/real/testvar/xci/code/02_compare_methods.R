rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

l <- read.csv('hca/real/xci/list.csv',header=F)

res <- sapply(setdiff(list.files('hca/real/testvar/result/EM_pm'),'perf'),function(celltype) {
  ## Lamian.pm
  
  lamiangene <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/',celltype,'/gender/gender_res.rds'))[[1]]
  lamiangene <- rownames(lamiangene)
  lamiansig <- rownames(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/EM_pm/',celltype,'/gender/differential_genes.csv'),as.is=T,row.names = 1))
  
  ## limma
  limma <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/limma/',celltype,'/gender_res.rds'))
  limmagene <- rownames(limma)[order(limma$P.Value)] ## all genes order
  limmasig <- rownames(limma)[limma$adj.P.Val < 0.05] ## significant genes
  ## tradeSeq
  tradeseq <- readRDS(paste0('hca/real/testvar/result/tradeSeq/',celltype,'/gender/testvar_res.rds'))
  tradeseqgene <- lapply(tradeseq,function(i) rownames(i)[order(i$P.Value,-i$waldStat)])
  tradeseqsig <- lapply(tradeseq,function(i) rownames(i)[i$adj.P.Val < 0.05])
  ## Lamian.chisq
  lamc <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/Lamian.chisq/',celltype,'/gender_res.rds'))[[1]]
  lamcgenes <- rownames(lamc)[order(lamc[,2])] ## all genes order
  lamcsig <- rownames(lamc)[lamc[,'fdr.chisq.overall'] < 0.05]   ## significant genes
  ## condiments
  res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/condiments/', celltype, '/gender/cond_gene_res.rds'))
  res[is.na(res[,3]),3] <- 1
  res$FDR <- p.adjust(res[,3],method='fdr')
  res = res[order(res[, 3],-abs(res[, 1])),]
  condiments_genes = rownames(res) ## all genes
  condiments_sig = rownames(res)[res[,'FDR']<0.05] ## significant genes
  ## monocle2 trajTest
  res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/monocle2_trajtest/', celltype, '/gender/res.rds'))
  res = res[order(res[, 1]),]
  monocle2trajtest_genes = rownames(res)
  monocle2trajtest_sig = rownames(res)[res[,'fdr'] < 0.05]
  ## monocle2 trajTest corrected by ourself (change null model)
  res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/monocle2_trajtest.corr/', celltype, '/gender/res.rds'))
  res = res[order(res[, 1]),]
  monocle2trajtestcorr_genes = rownames(res)
  monocle2trajtestcorr_sig = rownames(res)[res[,'fdr'] < 0.05]
  ## phenopath 100 epoch
  fit <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/phenopath100/', celltype, '/gender/fit_res.rds'))
  zscore <- abs(fit$m_beta[1,]/sqrt(fit$s_beta[1,])) 
  names(zscore) <- fit$feature_names
  pval <- pnorm(zscore,lower.tail = F)
  res = data.frame(score=zscore,pval=pval,fdr=p.adjust(pval,method='fdr'))
  res <- res[order(-res[,1],res[,2]),]
  phenopath_genes = rownames(res)
  phenopath_sig = rownames(res)[res[,'fdr'] < 0.05]
  
  ########
  genes <- list(Lamian=lamiangene)
  for (i in names(tradeseq))
    genes[[paste0('tradeSeq_',i)]]=tradeseqgene[[i]]
  genes[['limma']] <- limmagene
  genes[['condiments']] <- condiments_genes
  genes[['monocle2trajtest']] <- monocle2trajtest_genes
  genes[['phenopath']] <- phenopath_genes
  genes[['Lamian.chisq']] <- lamcgenes
  genes[['monocle2trajtest.corr']] <- monocle2trajtestcorr_genes
  
  sig <- list(Lamian=lamiansig)
  for (i in names(tradeseq))
    sig[[paste0('tradeSeq_',i)]]=tradeseqsig[[i]]
  sig[['limma']] <- limmagene
  sig[['condiments']] <- condiments_sig
  sig[['monocle2trajtest']] <- monocle2trajtest_sig
  sig[['phenopath']] <- phenopath_sig
  sig[['Lamian.chisq']] <- lamcsig
  sig[['monocle2trajtest.corr']] <- monocle2trajtestcorr_sig
  
  sapply(names(sig),function(n) {
    g <- sub(':.*','',sig[[n]])
    ng <- setdiff(sub(':.*','',genes[[n]]),g)
    t1 <- l[match(intersect(g,l[,1]),l[,1]),2]
    t2 <- l[match(intersect(ng,l[,1]),l[,1]),2]
    mean(t1 %in% c('E','Mostly E'))
  })
})

library(reshape2)
pd <- melt(res)
pdir <- 'hca/real/testvar/xci/plot/'
saveRDS(pd, paste0(pdir, 'compare_pd.rds'))

source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')

library(ggplot2)
pdf(paste0(pdir, '/compare.pdf'), height = 1.8, width = 3.6)
ggplot(pd,aes(x=Var1,y=value)) + geom_bar(stat='identity') + facet_wrap(~Var2) + coord_flip() + ylab('Percentage of differential genes escaping x-inactivation') + xlab('') + scale_y_continuous(breaks=c(0, 0.5, 1))
dev.off()

