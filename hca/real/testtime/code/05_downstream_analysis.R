library(here)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(grid)
setwd(here())
# path = 'lymph'
path = 'erythroid'
# path = 'monocyte'
rdir <- ddir <- paste0('hca/real/testtime/result/', path, '/')
pdir  <- paste0('hca/real/testtime/plot/', path, '/')
dir.create(pdir, recursive = T)
source('function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')

# ------------- #
# prepare data  #
# ------------- #
Res <- readRDS(paste0(ddir, '/testtime_res.rds'))
Res$populationFit <- readRDS(paste0(rdir, '/population_fit.rds'))
Res$cluster <- readRDS(paste0(rdir, '/cluster.rds'))

# ---------------------------- #
# downstream analysis pipeline #  
# ---------------------------- #

## --------------
## save diff gene
## --------------
allg <- names(Res$fdr[Res$fdr < 0.05])
res <- data.frame(fdr = Res$fdr[allg], fc = Res$foldchange[allg], pvalue = Res$pval[allg], stringsAsFactors = FALSE)
res <- res[order(res$fdr, abs(res$fc)), ]
write.csv(res, paste0(rdir, 'testtime_differential_genes.csv'))

## ---------------
## plotClusterMean
## ----------------
pdf(paste0(pdir, 'cluster_mean.pdf'), width = 5, height = 3.5)
plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'time')
dev.off()

## -----------
## GO analysis
## -----------
goRes <- GOEnrich(testobj = Res, type = 'time', version = 1, sep = ':.*')
saveRDS(goRes, paste0(rdir, 'goRes.rds'))

nn <- sapply(1:length(goRes), function(i){
  tmp <- goRes[[i]]
  # tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
  write.csv(tmp, paste0(rdir, 'cluster', i, '_GO.csv'))
  saveRDS(tmp, paste0(rdir, 'cluster', i, '_GO.rds'))
  print(str(tmp))
  return(0)
})

pdf(paste0(pdir, 'hm_GO_term.pdf'), width = 7.2, height = 3.5)
print(plotGOEnrich(goRes))
dev.off()


