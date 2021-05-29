library(here)
setwd(here())
ddir <- 'hca/real/build_from_tree_variability/result/monocyte/'
## load lamian, tradeSeq, limma results objects
lamian <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/su/res/Mod_Mi/numeric_res.rds')
lamsig <- rownames(lamian[[1]])[lamian[[1]][,1] < 0.05]
tradeseq <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/su/compres/tradeSeq/Mod_Mi/testvar_res.rds')
limma <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/su/compres/limma/Mod_Mi.rds')
limmasig <- rownames(limma)[limma$adj.P.Val < 0.05]
tradesig <- unique(unlist(sapply(tradeseq,function(i) rownames(i)[which(i$adj.P.Val < 0.05)])))
source('function/01_function.R')
lamian$populationFit <- getPopulationFit(lamian, rownames(lamian$expr), 'variable')

## load tradeSeq sce object and the counts
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
# setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
sce = readRDS('covid/Su_2020_Cell/testvar/useSaverImp/result/tradeSeq/Mod_Mi/sce.rds')
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
source('./function/01_function.R')
expr <- readRDS('covid/Su_2020_Cell/data/saver_log2norm_sub.rds')
expr <- expr[rowMeans(expr>0.1)>0.01, ]
ser <- readRDS('covid/Su_2020_Cell/data/CD8integrate.rds')
cnt <- as.matrix(ser@assays$RNA@counts)
cnt <- cnt[intersect(rownames(cnt),rownames(expr)), colnames(expr)]
cnt <- cnt[, colnames(sce)]

## plot example genes

### ============
### limma_excEM
### ============
targ <- intersect(rownames(lamian[[1]]),setdiff(limmasig,lamsig))
set.seed(12345)
selgene = sample(targ,100)
for (i in selgene) {
  print(i)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/limma_excEM/',i,'.pdf'))
  plotGeneSampleAndPopulation(lamian,i,variable='type',line.alpha=0.4)
  dev.off()
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/limma_excEM/', i, '_CI.pdf'), width = 1.9, height = 1.4)
  print(plotGeneCIAndPopulation(lamian, i, variable='type', ribbon.alpha=0.2, axis.text.blank = F, line.size = 0.9))
  dev.off()
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/limma_excEM/', i, '_samplefit.pdf'), width = 2, height = 1.5)
  plotGene(lamian, i, variable = 'type', ncol = 4, axis.text.blank = F, line.size = 0.2, continuous = T, use.palette = T)
  dev.off()
}

ag = c('SH2D3A', 'TMEM183A')
for (i in ag) {
  print(i)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/limma_excEM/',i,'.pdf'))
  plotGeneSampleAndPopulation(lamian,i,variable='type',line.alpha=0.4)
  dev.off()
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/limma_excEM/', i, '_CI.pdf'), width = 2.1, height = 1.4)
  print(plotGeneCIAndPopulation(lamian, i, variable='type', ribbon.alpha=0.2, axis.text.blank = F, line.size = 0.9))
  dev.off()
  
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/limma_excEM/', i, '_samplefit.pdf'), width = 2, height = 1.5)
  plotGene(lamian, i, variable = 'type', ncol = 4, axis.text.blank = F, line.size = 0.2, continuous = T, use.palette = T)
  dev.off()
}



### ===============
### tradeSeq_excEM
### ===============
targ <- intersect(rownames(lamian[[1]]),setdiff(tradesig,lamsig))
set.seed(12345)
selgene = sample(targ,100)
for (i in selgene) {
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/tradeseq_excEM/',i,'.pdf'),  width = 2, height = 1.5)
  plotGeneSampleAndPopulation(lamian,i,variable='type',line.alpha=0.4)
  dev.off()
  
  pdf(paste0(pdir, 'CI_', g, '.pdf'), width = 1.9, height = 1.4)
  print(plotGeneCIAndPopulation(lamian, g, variable='type', ribbon.alpha=0.2, axis.text.blank = F))
  dev.off()
  
  pdf(paste0(pdir, 'samplefit_', g, '.pdf'), width = 2, height = 1.5)
  plotGene(lamian, g, variable = 'type', ncol = 4, axis.text.blank = F, line.size = 0.2, continuous = T, use.palette = T)
  dev.off()
}

### ==================
### EM_extradeseqlimma
### ==================
targ <- setdiff(lamsig,c(tradesig,limmasig))
set.seed(12345)
selgene = sample(targ,100)
for (i in selgene) {
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/EM_extradeseqlimma/',i,'.pdf'),  width = 2, height = 1.5)
  plotGeneSampleAndPopulation(lamian,i,variable='type',line.alpha=0.4)
  dev.off()
}

## ============================
## plot selected tradeSeq_excEM
## ============================
ag = c('CCL4L2', 'CD300A')
for (g in ag) {
  pdf(paste0('covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/tradeseq_excEM/',g,'.pdf'), width = 2.9, height = 2)
  plot(plotSmoothers(sce, cnt, gene = g, size = 0.1, lwd = 1))
  dev.off()  
}
for (g in ag) {
  pdf(paste0('/Users/wenpinhou/Dropbox/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/plot/examplegene/tradeseq_excEM/',g,'_sample.pdf'), width = 2, height = 1.5)
  plotGene(lamian, g, variable = 'type', ncol = 4, axis.text.blank = T, line.size = 0.2, continuous = T, use.palette = T)
  dev.off()  
}

