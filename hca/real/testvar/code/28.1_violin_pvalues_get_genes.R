rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
# setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home/whou10/scratch16/whou10/')
# source('function/01_function.R')

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('resource/chrX_genename.rds')
u2 = readRDS('resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

for (celltype in setdiff(list.files('trajectory_variability/hca/real/testvar/result/EM_pm'),'perf')) {
  ## Lamian.pm
  lamian <- readRDS(paste0('trajectory_variability/hca/real/testvar/result/EM_pm/',celltype,'/gender/gender_res.rds'))[[1]]
  lamiangene <- rownames(lamian)[order(lamian[,2])] ## all genes order
  if (celltype=='monocyte') { ## significant genes
    lamiansig <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]   
  } else {
    lamiansig <- rownames(lamian)[lamian[,'pval.overall'] < 0.05]    
  }
  ## limma
  limma <- readRDS(paste0('trajectory_variability/hca/real/testvar/result/limma/',celltype,'/gender_res.rds'))
  limmagene <- rownames(limma)[order(limma$P.Value)] ## all genes order
  limmasig <- rownames(limma)[limma$adj.P.Val < 0.05] ## significant genes
  ## tradeSeq
  tradeseq <- readRDS(paste0('trajectory_variability/hca/real/testvar/result/tradeSeq/',celltype,'/gender/testvar_res.rds'))
  tradeseqgene <- lapply(tradeseq,function(i) rownames(i)[order(i$P.Value,-i$waldStat)])
  tradeseqsig <- lapply(tradeseq,function(i) rownames(i)[i$adj.P.Val < 0.05])
  ## Lamian.chisq
  lamc <- readRDS(paste0('trajectory_variability/hca/real/testvar/result/Lamian.chisq/',celltype,'/gender_res.rds'))[[1]]
  lamcgenes <- rownames(lamc)[order(lamc[,2])] ## all genes order
  lamcsig <- rownames(lamc)[lamc[,'fdr.chisq.overall'] < 0.05]   ## significant genes
  ## condiments
  res = readRDS(paste0('trajectory_variability/hca/real/testvar/result/condiments/', celltype, '/gender/cond_gene_res.rds'))
  res[is.na(res[,3]),3] <- 1
  res$FDR <- p.adjust(res[,3],method='fdr')
  res = res[order(res[, 3],-abs(res[, 1])),]
  condiments_genes = rownames(res) ## all genes
  condiments_sig = rownames(res)[res[,'FDR']<0.05] ## significant genes
  ## monocle2 trajTest
  res = readRDS(paste0('trajectory_variability/hca/real/testvar/result/monocle2_trajtest/', celltype, '/gender/res.rds'))
  res = res[order(res[, 1]),]
  monocle2trajtest_genes = rownames(res)
  monocle2trajtest_sig = rownames(res)[res[,'fdr'] < 0.05]
  ## monocle2 trajTest corrected by ourself (change null model)
  res = readRDS(paste0('trajectory_variability/hca/real/testvar/result/monocle2_trajtest.corr/', celltype, '/gender/res.rds'))
  res = res[order(res[, 1]),]
  monocle2trajtestcorr_genes = rownames(res)
  monocle2trajtestcorr_sig = rownames(res)[res[,'fdr'] < 0.05]
  ## phenopath 100 epoch
  fit <- readRDS(paste0('trajectory_variability/hca/real/testvar/result/phenopath100/', celltype, '/gender/fit_res.rds'))
  zscore <- abs(fit$m_beta[1,]/sqrt(fit$s_beta[1,])) 
  names(zscore) <- fit$feature_names
  pval <- pnorm(zscore,lower.tail = F)
  res = data.frame(score=zscore,pval=pval,fdr=p.adjust(pval,method='fdr'))
  res <- res[order(-res[,1],res[,2]),]
  phenopath_genes = rownames(res)
  phenopath_sig = rownames(res)[res[,'fdr'] < 0.05]
  
  
  #######
  ## save all significant genes
  siggene = list(
    Lamian.pm = lamiansig,
    Lamian.chisq = lamcsig,
    limma = limmasig,
    tradeSeq = tradeseqsig,
    condiments = condiments_sig,
    phenopath = phenopath_sig,
    monocle2trajTestCorr = monocle2trajtestcorr_sig,
    monocle2trajTest = monocle2trajtest_sig
  )
  saveRDS(siggene, paste0('trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_venn_significant_genes.rds'))
  
  ######## 2 (Lamian) + 3 (tradeseq patterns) + 6 (competing) = 11
  genes <- list(Lamian.pm=lamiangene)
  genes[['Lamian.chisq']] <- lamcgenes
  genes[['tradeSeqDT']] <- tradeseqgene[[1]]
  genes[['tradeSeqPT']] <- tradeseqgene[[2]]
  genes[['tradeSeqET']] <- tradeseqgene[[3]]
  genes[['tradeSeq']] <- unique(unlist(tradeseqgene))
  genes[['limma']] <- limmagene
  genes[['condiments']] <- condiments_genes
  genes[['monocle2trajTest']] <- monocle2trajtest_genes
  genes[['phenopath']] <- phenopath_genes
  genes[['monocle2trajTestCorr']] <- monocle2trajtestcorr_genes
  
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
    Lamian_excTradeSeq = setdiff(lamiangene, unique(unlist(tradeseqsig))), 
    Lamian_excCondiments = setdiff(lamiangene, condiments_sig),
    Lamian_excPhenopath = setdiff(lamiangene, phenopath_sig),
    Lamian_excMonocle2trajTest = setdiff(lamiangene, monocle2trajtest_sig),
    Lamian_excMonocle2trajTestCorr = setdiff(lamiangene, monocle2trajtestcorr_sig),
    Lamian_excTradeSeqDT = setdiff(lamiangene, unique(tradeseqsig[1])),
    Lamian_excTradeSeqPT = setdiff(lamiangene, unique(tradeseqsig[2])),
    Lamian_excTradeSeqET = setdiff(lamiangene, unique(tradeseqsig[3])),
    tradeSeqDT_excludeLamian = setdiff(tradeseqgene[[1]],lamiansig), 
    tradeSeqPT_excludeLamian = setdiff(tradeseqgene[[2]],lamiansig), 
    tradeSeqET_excludeLamian = setdiff(tradeseqgene[[3]],lamiansig), 
    tradeSeq_excLamian = setdiff(unique(c(unlist(tradeseqgene))), lamiansig),
    condiments_excLamian = setdiff(condiments_genes, lamiansig),
    phenopath_excLamian = setdiff(phenopath_genes, lamiansig),
    monocle2trajTest_excLamian = setdiff(monocle2trajtest_genes, lamiansig),
    monocle2trajTestCorr_excLamian = setdiff(monocle2trajtestcorr_genes, lamiansig))
  genes = c(genes, sig)
  saveRDS(genes, paste0('trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_violin_plotdata_genes.rds'))
}


