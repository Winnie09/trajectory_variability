rm(list=ls())
library(here)
library(ggplot2)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
source('function/01_function.R')

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
u2 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')
# u1 = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
# u2 = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')

Res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/result/EM_pm/monocyte/gender/gender_res.rds'))
lamian <- Res[[1]]
lamian <- rownames(lamian)[lamian[,'fdr.overall'] < 0.05]
limma <- readRDS(paste0('hca/real/testvar/result/limma/monocyte/gender_res.rds'))
limma <- rownames(limma)[limma$adj.P.Val < 0.05]
tradeseq <- readRDS(paste0('hca/real/testvar/result/tradeSeq/monocyte/gender/testvar_res.rds'))
tradeseq <- unique(unlist(sapply(tradeseq,function(i) {rownames(i)[i$adj.P.Val < 0.05]})))
tradeseq <- tradeseq[!is.na(tradeseq)]


a = setdiff(lamian, limma)
b = setdiff(a, tradeseq)
str(a)
str(b)


png('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/EM_pm/monocyte/gender/lamian_newgene.png', width = 600, height = 600)
plotGene(Res, b, plot.point = F, variable = 'gender', ncol = 4)
dev.off()


rn = rownames(Res$expr)
gene = c(rn[grepl('RCHY1', rn)], rn[grepl('LANCL1', rn)])
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/EM_pm/monocyte/gender/sampleFit_newgene3.pdf', width = 4.1, height = 2.1)
plotGene(Res, gene, plot.point = F, variable = 'gender', ncol = 2, continuous = F, sep = ':.*')
dev.off()

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/plot/EM_pm/monocyte/gender/sampleFit_newgene3_dots.pdf', width = 4.1, height = 2.1)
plotGene(Res, gene, plot.point = T, variable = 'gender', ncol = 2, continuous = F, line.alpha = 0, line.size = 0,  sep = ':.*', point.size = 0.5)
dev.off()
