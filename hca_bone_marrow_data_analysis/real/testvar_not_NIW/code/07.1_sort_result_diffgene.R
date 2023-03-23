rm(list=ls())
library(here)
library(ggplot2)
setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
# setwd(here())
source('function/01_function.R')
## read in gold standard 
agediff <- readRDS('agediff/result/age_diffgene.rds')
ChrX = readRDS('/Users/wenpinhou/Dropbox/resource/chrX_genename.rds')
ChrY = readRDS('/Users/wenpinhou/Dropbox/resource/chrY_genename.rds')
  
  
for (path in c('erythroid', 'lymph', 'monocyte')){
  ddir <- rdir <- paste0('hca/real/testvar/result/', path)
  pdir <- paste0('hca/real/testvar/plot/', path)
  dir.create(pdir)
  
  ## ---------------
  ## age difference
  ## ---------------
  Res <- readRDS(paste0(ddir, '/age/age_res.rds'))
  res = Res$statistics
  ## identify DDGType
  DDGType <- getDDGType(Res)
  res <- cbind(res, DDGType = DDGType[rownames(res)], TrurPositive = ifelse(sub(':.*','', rownames(res)) %in% agediff, 'TRUE', 'FALSE'))
  saveRDS(res, paste0(rdir, '/age/age_fdr_res.rds'))
  write.csv(res, paste0(rdir, '/age/age_fdr_res.csv'))
  sink(paste0(rdir, '/age/age_DDGType.txt'))
  table(DDGType)
  sink()
}
  ## plot DDG: all, meanSig, trendSig
  res <- res[order(res[,7], -res[,8]), ]
  png(paste0(pdir, '/age/age_diffgene.png'), width = 1700, height = 1500, res = 200)
  plotGenePopulation(Res, rownames(res)[1:min(25, nrow(res))], type = 'variable', sep = ':.*')
  dev.off()

  res <- res[order(res[,1], -res[,2]), ]
  int <- intersect(rownames(res)[res[,1]<0.05], names(DDGType)[DDGType == 'meanSig'])
  png(paste0(pdir, '/age/age_diffgene_meanSig.png'), width = 1700, height = 1500, res = 200)
  plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
  dev.off()

  res <- res[order(res[,4], -res[,5]), ]
  int <- intersect(rownames(res)[res[,4]<0.05], names(DDGType)[DDGType == 'trendSig'])
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
  DDGType <- getDDGType(Res)
  res <- cbind(res, DDGType = DDGType[rownames(res)], 
               inChrX = ifelse(sub(':.*', '', rownames(res)) %in% ChrX, 'TRUE', 'FALSE'),
               inChrY = ifelse(sub(':.*', '', rownames(res)) %in% ChrY, 'TRUE', 'FALSE'))
  saveRDS(res, paste0(rdir, '/gender/gender_fdr_res.rds'))
  write.csv(res, paste0(rdir, '/gender/gender_fdr_res.csv'))

  sink(paste0(rdir, '/gender/gender_DDGType.txt'))
  table(DDGType)
  sink()
  
  res <- res[order(res[,7], -res[,8]), ]
  png(paste0(pdir, '/gender/gender_diffgene.png'), width = 1700, height = 1500, res = 200)
  plotGenePopulation(Res, rownames(res)[1:min(25,nrow(res))], type = 'variable', sep = ':.*')
  dev.off()
  
  res <- res[order(res[,1], -res[,2]), ]
  int <- intersect(rownames(res)[res[,1]<0.05], names(DDGType)[DDGType == 'meanSig'])
  if (length(int) > 0){
    png(paste0(pdir, '/gender/gender_diffgene_meanSig.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
    dev.off()
  }
    
  res <- res[order(res[,4], -res[,5]), ]
  int <- intersect(rownames(res)[res[,4]<0.05], names(DDGType)[DDGType == 'trendSig'])
  if (length(int) > 0){
    png(paste0(pdir, '/gender/gender_diffgene_trendSig.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
    dev.off()
  }
}

