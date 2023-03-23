rm(list=ls())
library(here)
library(ggplot2)
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability')
setwd(here())
source('function/01_function.R')
## read in gold standard 
agediff <- readRDS('agediff/result/age_diffgene.rds')
ChrX = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrX_genename.rds')
ChrY = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/chrY_genename.rds')

for (path in c('erythroid', 'lymph', 'monocyte')){
  ddir <- rdir <- paste0('hca/real/testvar/result/EM_pm/', path)
  pdir <- paste0('hca/real/testvar/plot/EM_pm/', path)
  dir.create(paste0(pdir, '/age'), recursive = T)
  dir.create(paste0(pdir, '/gender'), recursive = T)
  ## ---------------
  ## age difference
  ## ---------------
  Res <- readRDS(paste0(ddir, '/age/age_res.rds'))
  res = Res$statistics
  ## identify DDGType
  DDGType <- getDDGType(Res)
  sink(paste0(rdir, '/age/age_DDGType.txt'))
  table(DDGType)
  sink()
  
  ## plot DDG
  if (sum(res[,1] < 0.05) > 0){
    res <- cbind(res, DDGType = DDGType[rownames(res)], TrurPositive = ifelse(sub(':.*','', rownames(res)) %in% agediff, 'TRUE', 'FALSE'))
    saveRDS(res, paste0(rdir, '/age/age_fdr_res.rds'))
    write.csv(res, paste0(rdir, '/age/age_fdr_res.csv'))
    ## plot DDG in different types
    res <- Res$statistics
    res <- res[res[,1]<0.05, ,drop=F]
    res <- res[order(res[,'fdr.overall'], -res[,'z.overall']), ]
    DDGType <- DDGType[rownames(res)]
    png(paste0(pdir, '/age/age_diffgene.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, rownames(res)[1:min(25, nrow(res))], type = 'variable', sep = ':.*')
    dev.off()
    
    for (i in setdiff(unique(DDGType), 'nonDDG')){
      int <- names(DDGType)[DDGType == i]
      if (length(int) > 0){
        png(paste0(pdir, '/gender/gender_diffgene_', i, '.png'), width = 1700, height = 1500, res = 200)
        plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
        dev.off()
      }  
    }
    
  }
  ### ------------------
  ### gender difference
  ### ------------------
  Res <- readRDS(paste0(ddir, '/gender/gender_res.rds'))
  res <- Res$statistics
  ## identify DDGType
  DDGType <- getDDGType(Res)
  sink(paste0(rdir, '/gender/gender_DDGType.txt'))
  table(DDGType)
  sink()
  ## plot DDG
  if (sum(res[,1] < 0.05) > 0){
    res <- cbind(res, DDGType = DDGType[rownames(res)], 
                 inChrX = ifelse(sub(':.*', '', rownames(res)) %in% ChrX, 'TRUE', 'FALSE'),
                 inChrY = ifelse(sub(':.*', '', rownames(res)) %in% ChrY, 'TRUE', 'FALSE'))
    saveRDS(res, paste0(rdir, '/gender/gender_fdr_res.rds'))
    write.csv(res, paste0(rdir, '/gender/gender_fdr_res.csv'))
    ## plot DDG in different types
    res = Res$statistics
    res <- res[res[,1] < 0.05, ,drop=F]
    res <- res[order(res[,'fdr.overall'], -res[,'z.overall']), ]
    DDGType <- DDGType[rownames(res)]
    png(paste0(pdir, '/gender/gender_diffgene.png'), width = 1700, height = 1500, res = 200)
    plotGenePopulation(Res, rownames(res)[1:min(25,nrow(res))], type = 'variable', sep = ':.*')
    dev.off()
    
    for (i in setdiff(unique(DDGType), 'nonDDG')){
      int <- names(DDGType)[DDGType == i]
      if (length(int) > 0){
        png(paste0(pdir, '/gender/gender_diffgene_', i, '.png'), width = 1700, height = 1500, res = 200)
        plotGenePopulation(Res, int[1:min(25, length(int))], type = 'variable', sep = ':.*')
        dev.off()
      }  
    }
  } 
}

