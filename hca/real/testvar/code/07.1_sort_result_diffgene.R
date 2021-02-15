rm(list=ls())
library(here)
library(ggplot2)
setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# setwd(here())
source('function/01_function.R')

for (path in c('erythroid', 'lymph', 'monocyte')){
  ddir <- rdir <- paste0('hca/real/testvar/result/', path)
  pdir <- paste0('hca/real/testvar/plot/', path)
  dir.create(pdir)
  
  ## ---------------
  ## age difference
  ## ---------------
  Res <- readRDS(paste0(ddir, '/age/age_res.rds'))
  res = Res$statistics
  ## identify DEGType 
  DEGType <- getDEGType(Res)
  res <- cbind(res, DEGType = DEGType[rownames(res)])
  saveRDS(res, paste0(rdir, '/age_fdr_res.rds'))
  write.csv(res, paste0(rdir, '/age_fdr_res.csv'))
  sink(paste0(rdir, '/age/age_DEGType.txt'))
  table(DEGType)
  sink()
  
  ## plot DEG: all, meanSig, trendSig
  res <- res[order(res[,7], -res[,8]), ]
  png(paste0(pdir, '/age/age_diffgene.png'), width = 1700, height = 1500, res = 200)
  plotGenePopulation(Res, rownames(res)[1:min(25, nrow(res))], type = 'variable', sep = ':.*')
  dev.off()
  
  res <- res[order(res[,1], -res[,2]), ]
  int <- intersect(rownames(res)[res[,1]<0.05], names(DEGType)[DEGType == 'meanSig'])
  png(paste0(pdir, '/age/age_diffgene_meanSig.png'), width = 1700, height = 1500, res = 200)
  plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
  dev.off()

  res <- res[order(res[,4], -res[,5]), ]
  int <- intersect(rownames(res)[res[,4]<0.05], names(DEGType)[DEGType == 'trendSig'])
  if (length(int) > 0){
    png(paste0(pdir, '/age/age_diffgene_trendSig.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
    dev.off()
  }
    
  ### ------------------
  ### gender difference
  ### ------------------
  Res <- readRDS(paste0(ddir, '/gender/gender_res.rds'))
  res <- Res$statistics
  DEGType <- getDEGType(Res)
  res <- cbind(res, DEGType = DEGType[rownames(res)])
  saveRDS(res, paste0(rdir, '/gender/gender_fdr_res.rds'))
  write.csv(res, paste0(rdir, '/gender/gender_fdr_res.csv'))
  
  sink(paste0(rdir, '/gender_DEGType.txt'))
  table(DEGType)
  sink()
  
  res <- res[order(res[,7], -res[,8]), ]
  png(paste0(pdir, '/gender/gender_diffgene.png'), width = 1700, height = 1500, res = 200)
  plotGenePopulation(Res, rownames(res)[1:min(25,nrow(res))], type = 'variable', sep = ':.*')
  dev.off()
  
  res <- res[order(res[,1], -res[,2]), ]
  int <- intersect(rownames(res)[res[,1]<0.05], names(DEGType)[DEGType == 'meanSig'])
  if (length(int) > 0){
    png(paste0(pdir, '/gender/gender_diffgene_meanSig.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
    dev.off()
  }
    
  res <- res[order(res[,4], -res[,5]), ]
  int <- intersect(rownames(res)[res[,4]<0.05], names(DEGType)[DEGType == 'trendSig'])
  if (length(int) > 0){
    png(paste0(pdir, '/gender/gender_diffgene_trendSig.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
    dev.off()
  }
}

