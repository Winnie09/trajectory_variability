rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# source('/Users/wenpinhou/Dropbox/resource/myfunc/01_function.R')
source('function/01_function.R')
for (dataType in seq(1, 4)) {
  print(dataType)
  ddir <- 'hca/simu/testvar/addMultiSignalUsingExpr/'
  rdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/cluster/',dataType, '/')
  pdir <- paste0('hca/simu/testvar/addMultiSignalUsingExpr/plot/', dataType, '/')
  dir.create(pdir)
  dir.create(rdir, recursive = T)
  Res <- readRDS(paste0('hca/simu/testvar/addMultiSignalUsingExpr/result/EM_NOT_centered/',dataType,'.rds'))
  res <- Res$statistics 
  ## diffType
  library(limma)
  diffType <- cbind(meanDiff = (res$meanDiff.fdr < 0.05),
                    trendDiff = (res$trendDiff.fdr < 0.05),
                    both = (res$both.fdr < 0.05),
                    meanTrue = (rownames(res) %in% meangene),
                    trendTrue = (rownames(res) %in% trendgene),
                    meantrendTrue = (rownames(res) %in% meantrendgene))
  pdf(paste0(pdir, 'venn.pdf'), width = 15, height = 3.5)
  par(mfrow=c(1,4))
  vennDiagram(diffType[,1:3])
  vennDiagram(diffType[,c(1,4:6)])
  vennDiagram(diffType[,c(2,4:6)])
  vennDiagram(diffType[,c(3,4:6)])
  dev.off()
  
  res = Res$statistics
  diffType <- cbind(meanOnly = ((res$meanDiff.fdr < 0.05) & (res$both.fdr < 0.05) & (res$trendDiff.fdr > 0.05) | ((res$both.fdr > 0.05) & (res$meanDiff.fdr < 0.05) & (res$trendDiff.fdr > 0.05))),
                    trendOnly = (((res$trendDiff.fdr < 0.05) & (res$both.fdr < 0.05) & (res$meanDiff.fdr > 0.05))| ((res$both.fdr > 0.05) & (res$trendDiff.fdr < 0.05) & (res$meanDiff.fdr > 0.05))),
                    both = ((res$both.fdr < 0.05) & (res$meanDiff.fdr < 0.05) & (res$trendDiff.fdr < 0.05)),
                    meanTrue = (rownames(res) %in% meangene),
                    trendTrue = (rownames(res) %in% trendgene),
                    meantrendTrue = (rownames(res) %in% meantrendgene))
  pdf(paste0(pdir, 'venn_only.pdf'), width = 15, height = 3.5)
  par(mfrow=c(1,4))
  vennDiagram(diffType[,1:3])
  vennDiagram(diffType[,c(1,4:6)])
  vennDiagram(diffType[,c(2,4:6)])
  vennDiagram(diffType[,c(3,4:6)])
  dev.off()
  
  ## classify to four types
  diffType <- cbind(meanOnly = ((res$meanDiff.fdr < 0.05) & (res$both.fdr < 0.05) & (res$trendDiff.fdr > 0.05)),
                    trendOnly = ((res$trendDiff.fdr < 0.05) & (res$both.fdr < 0.05) & (res$meanDiff.fdr > 0.05)),
                    both = ((res$both.fdr < 0.05) & (res$meanDiff.fdr < 0.05) & (res$trendDiff.fdr < 0.05)),
                    unknown = ((res$both.fdr < 0.05) & (res$meanDiff.fdr > 0.05) & (res$trendDiff.fdr > 0.05)),
                    meanTrue = (rownames(res) %in% meangene),
                    trendTrue = (rownames(res) %in% trendgene),
                    meantrendTrue = (rownames(res) %in% meantrendgene))
  rownames(diffType) = rownames(res)
  pdf(paste0(pdir, 'venn_only_withinDiff.pdf'), width = 19, height = 3.5)
  par(mfrow=c(1,5))
  vennDiagram(diffType[,1:4])
  vennDiagram(diffType[,c(1,5:7)])
  vennDiagram(diffType[,c(2,5:7)])
  vennDiagram(diffType[,c(3,5:7)])
  vennDiagram(diffType[,c(4,5:7)])
  dev.off()
  pdf(paste0(pdir, 'meanOnlyT_meanTrueT.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & diffType[,5]][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'meanOnlyT_meanTrueF.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & (!diffType[,5])][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'meanOnlyT_trendTrueT.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & diffType[,6]][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'meanOnlyT_trendTrueF.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & (!diffType[,6])][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
   pdf(paste0(pdir, 'meanOnlyT_meantrendTrueT.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & diffType[,7]][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
   pdf(paste0(pdir, 'meanOnlyT_meantrendTrueF.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & (!diffType[,7])][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  diffType <- cbind(diff = (res$both.fdr < 0.05),
                    meanDiff = ((res$meanDiff.fdr < 0.05) & (res$both.fdr < 0.05) ),
                    trendDiff = ((res$trendDiff.fdr < 0.05) & (res$both.fdr < 0.05)),
                    meanTrue = (rownames(res) %in% meangene),
                    trendTrue = (rownames(res) %in% trendgene),
                    meantrendTrue = (rownames(res) %in% meantrendgene))
  rownames(diffType) = rownames(res)
  pdf(paste0(pdir, 'venn_withinDiff_either.pdf'), width = 15, height = 3.5)
  par(mfrow=c(1,4))
  vennDiagram(diffType[,1:3])
  vennDiagram(diffType[,c(1,4:6)])
  vennDiagram(diffType[,c(2,4:6)])
  vennDiagram(diffType[,c(3,4:6)])
  dev.off()
  
  pdf(paste0(pdir, 'bothT_meanT_trendF.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & diffType[,2] & (!diffType[,3])][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'bothT_meanF_trendT.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & (!diffType[,2]) & diffType[,3]][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'bothT_meanF_trendF.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & (!diffType[,2]) & (!diffType[,3])][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'bothT_meanT_trendT.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[diffType[,1] & diffType[,2] & diffType[,3]][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'bothF_meanT_trendF.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[(!diffType[,1]) & diffType[,2] & (!diffType[,3])][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  pdf(paste0(pdir, 'bothF_meanT_trendT.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[(!diffType[,1]) & diffType[,2] & diffType[,3]][1:9], variable = 'group', sep = ':.*')
  dev.off()
  
  # pdf(paste0(pdir, 'bothF_meanF_trendT.pdf'), width = 7, height = 7) ## NA
  # plotGene(Res, gene = rownames(diffType)[(!diffType[,1]) & (!diffType[,2]) & diffType[,3]][1:9], variable = 'group', sep = ':.*')
  # dev.off()
  
  pdf(paste0(pdir, 'bothF_meanF_trendF.pdf'), width = 7, height = 7)
  plotGene(Res, gene = rownames(diffType)[(!diffType[,1]) & (!diffType[,2]) & (!diffType[,3])][1:9], variable = 'group', sep = ':.*')
  dev.off()
}
