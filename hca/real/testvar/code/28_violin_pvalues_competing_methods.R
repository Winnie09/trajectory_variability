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

for (celltype in setdiff(list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm'),'perf')) {
  ## Lamian.pm
  lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/',celltype,'/gender/gender_res.rds'))[[1]]
  lamiangene <- rownames(lamian)[order(lamian[,2])] ## all genes order
  if (celltype=='monocyte') { ## significant genes
    lamiansig <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]   
  } else {
    lamiansig <- rownames(lamian)[lamian[,'pval.overall'] < 0.05]    
  }
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
    Lamian_excTradeseq = setdiff(lamiangene, unique(unlist(tradeseqgene))), ## no genes
    Lamian_excCondiments = setdiff(lamiangene, condiments_sig),
    Lamian_excPhenopath = setdiff(lamiangene, phenopath_sig),
    Lamian_excMonocle2 = setdiff(lamiangene, monocle2trajtest_sig),
    tradeSeqDT_excludeLamian = setdiff(tradeseqgene[[1]],lamiansig), 
    tradeSeqPT_excludeLamian = setdiff(tradeseqgene[[2]],lamiansig), 
    tradeSeqET_excludeLamian = setdiff(tradeseqgene[[2]],lamiansig), 
    condiments_excLamian = setdiff(condiments_genes, lamiansig),
    phenopath_excLamian = setdiff(phenopath_genes, lamiansig),
    monocle2_excLamian = setdiff(monocle2trajtest_genes, lamiansig),
    monocle2corr_excLamian = setdiff(monocle2trajtestcorr_genes, lamiansig))
  genes = c(genes, sig)
  saveRDS(genes, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_violin_plotdata_genes.rds'))
  
  allg <- sub(':.*','',rownames(tradeseq[[1]]))
  
  permud <- reald <- NULL
  # for (met in names(sig)) {
  #   gn <- sig[[met]]
  #   gn <- sub(':.*','',gn)
  #   v1 <- mean(gn %in% u1)
  #   v2 <- mean(gn %in% u2)
  #   
  #   ##### permute reported gene order
  #   v1_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
  #     set.seed(myseed+100)
  #     gnpm = sample(allg,length(gn))
  #     mean(gnpm %in% u1)
  #   },mc.cores=detectCores()))
  #   
  #   v2_pm <- unlist(mclapply(seq(1,1e4), function(myseed){
  #     set.seed(myseed+100)
  #     gnpm = sample(allg,length(gn))
  #     mean(gnpm %in% u2)
  #   },mc.cores=detectCores()))
  #   
  #   permud <- rbind(permud,data.frame(per = c(v1_pm,v2_pm), type=rep(c('chrX','chrY'),each=1e4),method=met, stringsAsFactors = FALSE))
  #   reald <- rbind(reald,data.frame(per=c(v1,v2),type=c('chrX','chrY'),pvalue=c(mean(v1_pm >= v1),mean(v2_pm >= v2)),method=met,stringsAsFactors = F))
  # }
  
  for (met in names(genes)) {
    gn <- genes[[met]]
    allg <- sub(':.*','',gn)
    v1 <- cumsum(allg %in% u1)/seq(1,length(allg))
    v1 <- mean(v1[seq(1,length(v1)) %% 10 == 0])
    v2 <- cumsum(allg %in% u2)/seq(1,length(allg))
    v2 <- mean(v2[seq(1,length(v2)) %% 10 == 0])
    # saveRDS(v1, paste0(rdir, path, '/gender/gender_chrX_overlap.rds'))
    # saveRDS(v2, paste0(rdir, path, '/gender/gender_chrY_overlap.rds'))
    
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
  saveRDS(list(permud = permud, reald = reald), paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin_plotdata.rds'))
  
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin_all.pdf'),width=13,height=4)
  print(ggplot() + 
          geom_violin(data=permud,aes(x=method,y=per,col=type)) + 
          geom_point(data=reald,aes(x=method,y=per,col=type),size=1) + 
          geom_text(data=reald,aes(x=method,y=max(reald$per)*1.3,label=pvalue), size = 10*5/14)+
          theme_compact() + facet_wrap(~type) + 
          coord_flip(ylim=c(0,max(reald$per)*1.5)) + 
          xlab('') + 
          ylab('Proportion') + 
          scale_color_manual(values=c('chrX'=brewer.pal(3,'Pastel1')[1],'chrY'=brewer.pal(3,'Pastel1')[2])) +
          theme(legend.position = 'bottom',strip.background = element_blank(),strip.text = element_text(size=10),legend.title = element_blank(), text = element_text(size = 10)) + 
          scale_y_continuous(breaks=c(0,round(max(reald$per)*0.6,2),round(max(reald$per)*1.2,2)))
  )
  dev.off()
  
  
  # ## _pvalue_violin_all_exc3.pdf
  # excComp <- c('Lamian_excPhenopath', 'Lamian_excMonocle2', 'Lamian_excTradeseq')
  # 
  # ## _pvalue_violin_all_exc4.pdf
  # excComp <- c('Lamian_excPhenopath', 'Lamian_excMonocle2', 'Lamian_excTradeseq', 'Lamian_excCondiments')
  # 
  # ## _pvalue_violin_all_exc5.pdf
  # excComp <- c('Lamian_excPhenopath', 'Lamian_excMonocle2', 'Lamian_excTradeseq', 'Lamian_excCondiments', 'Lamian_excLimma')
  
  excComp <- c('Lamian_excPhenopath', 'Lamian_excMonocle2', 'Lamian_excTradeseq', 'Lamian_excCondiments', 'monocle2_excLamian', 'phenopath_excLamian')
  permud <- permud[!permud[,3] %in% excComp, ]
  reald <- reald[!reald[,4] %in% excComp, ]
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/perf/',celltype,'_pvalue_violin_all_exc6.pdf'),width=6.5,height=3.8)
  print(ggplot() + 
          geom_violin(data=permud,aes(x=method,y=per,col=type)) + 
          geom_point(data=reald,aes(x=method,y=per,col=type),size=1) + 
          geom_text(data=reald[reald[,2] == 'chrX', ],aes(x=method,y=max(reald$per)*1.3,label=pvalue), size = 10*5/14)+
          geom_text(data=reald[reald[,2] == 'chrY', ],aes(x=method,y=max(reald$per)*1.3 + 0.007,label=pvalue), size = 10*5/14)+
          theme_compact() + #facet_wrap(~type) + 
          coord_flip(ylim=c(0,max(reald$per)*1.5)) + 
          xlab('') + 
          ylab('Proportion') + 
          scale_color_manual(values=c('chrX'=brewer.pal(3,'Pastel1')[1],'chrY'=brewer.pal(3,'Pastel1')[2])) +
          theme(legend.position = 'bottom',strip.background = element_blank(),strip.text = element_text(size=7),legend.title = element_blank(), text = element_text(size = 7)) + 
          scale_y_continuous(breaks=c(0,round(max(reald$per)*0.6,2),round(max(reald$per)*1.2,2)), limits = c(0, max(reald$per) + 0.1))
  )
  dev.off()
}




