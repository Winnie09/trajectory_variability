library(ggplot2)
l <- list()
lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/Mod_Mi/numeric_res.rds'))[[1]]
l[['lamian.pm']] <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]
gn <- rownames(lamian)

lamian <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/Lamian.chisq/Mod_Mi/numeric_res.rds'))[[1]]
l[['lamian.chisq']] <- rownames(lamian)[lamian[,'fdr.chisq.overall'] < 0.05]


limma <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/limma/Mod_Mi.rds')
l[['limma']] <- rownames(limma)[limma$adj.P.Val < 0.05]
tradeseq <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/tradeSeq/Mod_Mi/testvar_res.rds')
for (n in names(tradeseq)) tradeseq[[n]]$adj.P.Val[is.na(tradeseq[[n]]$adj.P.Val)] <- 1
for (n in names(tradeseq)) l[[paste0('tradeseq_',n)]] <- rownames(tradeseq[[n]])[tradeseq[[n]]$adj.P.Val<0.05]

condiments <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/condiments/Mod_Mi/cond_gene_res.rds')
l[['condiments']] <- rownames(condiments)[p.adjust(condiments[,3],method='fdr') < 0.05]
monocle2_trajtest <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/monocle2_trajtest/Mod_Mi/res.rds')
l[['monocle2_trajtest']] <- rownames(monocle2_trajtest)[monocle2_trajtest$fdr < 0.05]
monocle2_trajtest.corr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/monocle2_trajtest.corr/Mod_Mi/res.rds')
l[['monocle2_trajtest.corr']] <- rownames(monocle2_trajtest.corr)[monocle2_trajtest.corr$fdr < 0.05]

phenopath <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/phenopath3/Mod_Mi/sig_res.rds')
tmp <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/phenopath3/Mod_Mi/fit_res.rds')
l[['phenopath']] <- tmp$feature_names[phenopath[,1]]

m <- matrix(0,nrow=length(gn),ncol=length(l),dimnames = list(gn,names(l)))
for (n in names(l)) m[l[[n]],n] <- 1


colSums(m)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/perf/pd_heatmap.rds')


pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/perf/heatmap.pdf',width=10,height=15)
heatmap(m,scale = 'none',labRow = NA)
dev.off()
