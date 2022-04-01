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

## Lamian.pm
lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex/pc2/lamian_pm_allcores.rds'))[[1]]
lamiangene <- rownames(lamian)[order(lamian[,2])] ## all genes order
lamiansig <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]   ## significant genes

## limma
limma <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex_noBatch/pc2/limma.rds'))
limmagene <- rownames(limma)[order(limma$P.Value)] ## all genes order
limmasig <- rownames(limma)[limma$adj.P.Val < 0.05] ## significant genes
## tradeSeq
# tradeseq <- readRDS(paste0('hca/real/testvar/result/tradeSeq/',celltype,'/gender/testvar_res.rds'))
# tradeseqgene <- lapply(tradeseq,function(i) rownames(i)[order(i$P.Value,-i$waldStat)])
# tradeseqsig <- lapply(tradeseq,function(i) rownames(i)[i$adj.P.Val < 0.05])
## Lamian.chisq
lamc <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex/pc2/lamian_chisq.rds'))[[1]]
lamcgenes <- rownames(lamc)[order(lamc[,2])] ## all genes order
lamcsig <- rownames(lamc)[lamc[,'fdr.chisq.overall'] < 0.05]   ## significant genes
## condiments
# res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/condiments/', celltype, '/gender/cond_gene_res.rds'))
# res[is.na(res[,3]),3] <- 1
# res$FDR <- p.adjust(res[,3],method='fdr')
# res = res[order(res[, 3],-abs(res[, 1])),]
# condiments_genes = rownames(res) ## all genes
# condiments_sig = rownames(res)[res[,'FDR']<0.05] ## significant genes
## monocle2 trajTest
res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex_noBatch/pc2/monocle2_trajtest.rds'))
res = res[order(res[, 1]),]
monocle2trajtest_genes = rownames(res)
monocle2trajtest_sig = rownames(res)[res[,'fdr'] < 0.05]
## monocle2 trajTest corrected by ourself (change null model)
res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/res/sex_noBatch/pc2/monocle2_trajtest_corr.rds'))
res = res[order(res[, 1]),]
monocle2trajtestcorr_genes = rownames(res)
monocle2trajtestcorr_sig = rownames(res)[res[,'fdr'] < 0.05]
## phenopath 100 epoch
# fit <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/phenopath100/', celltype, '/gender/fit_res.rds'))
# zscore <- abs(fit$m_beta[1,]/sqrt(fit$s_beta[1,])) 
# names(zscore) <- fit$feature_names
# pval <- pnorm(zscore,lower.tail = F)
# res = data.frame(score=zscore,pval=pval,fdr=p.adjust(pval,method='fdr'))
# res <- res[order(-res[,1],res[,2]),]
# phenopath_genes = rownames(res)
# phenopath_sig = rownames(res)[res[,'fdr'] < 0.05]

########
genes <- list(Lamian=lamiangene)
# for (i in names(tradeseq))
#   genes[[paste0('tradeSeq_',i)]]=tradeseqgene[[i]]
genes[['limma']] <- limmagene
#genes[['condiments']] <- condiments_genes
genes[['monocle2trajtest']] <- monocle2trajtest_genes
genes[['monocle2trajtestcorr']] <- monocle2trajtestcorr_genes
#genes[['phenopath']] <- phenopath_genes
genes[['Lamian.chisq']] <- lamcgenes
genes[['monocle2trajtest.corr']] <- monocle2trajtestcorr_genes

## ================================== 
## old way of calculating overlap score
## ================================== 
# sig <- list(
#   Lamian_excLimma = setdiff(lamiansig,limmasig),
#   Lamian_excTradeseq = setdiff(lamiansig, unique(unlist(tradeseqgene))), ## no genes
#   Lamian_excCondiments = setdiff(lamiansig, condiments_sig),
#   Lamian_excPhenopath = setdiff(lamiansig, phenopath_sig),
#   Lamian_excMonocle2 = setdiff(lamiansig, monocle2trajtest_sig),
#   tradeSeq_excludeLamian = setdiff(unlist(tradeseqsig),lamiansig), 
#   condiments_excLamian = setdiff(condiments_sig, lamiansig),
#   phenopath_excLamian = setdiff(phenopath_sig, lamiansig),
#   monocle2_excLamian = setdiff(monocle2trajtest_sig, lamiansig))
## ==============================================
## new way of calculating overlap score: 20211101
## ==============================================
sig <- list(
  Lamian_excLimma = setdiff(lamiangene,limmasig),
  Lamian_excMonocle2 = setdiff(lamiangene, monocle2trajtest_sig),
  Lamian_excMonocle2corr = setdiff(lamiangene, monocle2trajtestcorr_sig),
  monocle2_excLamian = setdiff(monocle2trajtest_genes, lamiansig),
  monocle2corr_excLamian = setdiff(monocle2trajtestcorr_genes, lamiansig))
genes = c(genes, sig)
saveRDS(genes, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/perf/violin_plotdata_genes.rds'))

allg <- lamiangene

permud <- reald <- NULL

for (met in names(genes)) {
  gn <- genes[[met]]
  allg <- sub(':.*','',gn)
  v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
  v1 <- mean(v1[seq(1,length(v1)) %% 10 == 0])
  v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
  v2 <- mean(v2[seq(1,length(v2)) %% 10 == 0])
  ##### permute reported gene order
  v1_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed+100)
    allg.pm = sample(allg)
    tmp = cumsum(allg.pm %in% u1)/seq(1,length(allg.pm))
    mean(tmp[seq(1,length(tmp)) %% 10 == 0])
  },mc.cores=detectCores()))
  
  v2_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
    set.seed(myseed+100)
    allg.pm = sample(allg)
    tmp = cumsum(allg.pm %in% u2)/seq(1,length(allg.pm))
    mean(tmp[seq(1,length(tmp)) %% 10 == 0])
  },mc.cores=detectCores()))
  
  permud <- rbind(permud,data.frame(per = c(v1_pm,v2_pm), type=rep(c('chrX','chrY'),each=1e4),method=met, stringsAsFactors = FALSE))
  reald <- rbind(reald,data.frame(per=c(v1,v2),type=c('chrX','chrY'),pvalue=c(mean(v1_pm >= v1),mean(v2_pm >= v2)),method=met,stringsAsFactors = F))
}
saveRDS(list(permud = permud, reald = reald), paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/perf/pvalue_violin_plotdata.rds'))

source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/perf/pvalue_violin_all.pdf',width=7.2,height=3.6)
print(ggplot() + 
          geom_violin(data=permud,aes(x=method,y=per,col=type), scale = 'width') + 
          geom_point(data=reald,aes(x=method,y=per,col=type),size=1) + 
          geom_text(data=reald[reald[,2] == 'chrX', ],aes(x=method,y=max(reald$per)*1.3,label=pvalue), size = 10*5/14)+
          geom_text(data=reald[reald[,2] == 'chrY', ],aes(x=method,y=max(reald$per)*1.3 + 0.02,label=pvalue), size = 10*5/14)+
          theme_compact() + #facet_wrap(~type) + 
          coord_flip(ylim=c(0,max(reald$per)*1.5)) + 
          xlab('') + 
          ylab('Proportion') + 
          scale_color_manual(values=c('chrX'=brewer.pal(3,'Pastel1')[1],'chrY'=brewer.pal(3,'Pastel1')[2])) +
          theme(legend.position = 'bottom',strip.background = element_blank(),strip.text = element_text(size=7),legend.title = element_blank(), text = element_text(size = 7)) + 
          scale_y_continuous(breaks=c(0,round(max(reald$per)*0.6,2),round(max(reald$per)*1.2,2)), limits = c(0, max(reald$per) + 0.1))
  )
dev.off()

