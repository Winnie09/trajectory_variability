rm(list=ls())
library(here)
# setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home/whou10/scratch16/whou10/')
# source('function/01_function.R')

## read in gold standard Sex difference genes (chrX, chrY)
u1 = readRDS('resource/chrX_genename.rds')
u2 = readRDS('resource/chrY_genename.rds')
a = readRDS('monocyte_venn_significant_genes.rds')
s = setdiff(a[[1]], a[[3]])
s = sub(':.*', '', s)

tradeseq <- readRDS(paste0('trajectory_variability/hca/real/testvar/result/tradeSeq/monocyte/gender/testvar_res.rds'))
allg <- sub(':.*','',rownames(tradeseq[[1]]))

chrx <- allg %in% u1 + 0
report <- allg %in% s + 0
tab = table(chrx, report)
tab
fisher.test(tab)

#     report
# chrx     0     1
#    0 11755    34
#    1   345     4
# 	Fisher's Exact Test for Count Data
# 
# data:  tab
# p-value = 0.02298
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.027426 11.327795
# sample estimates:
# odds ratio 
#   4.007672   

chry <- allg %in% u2 + 0

tab = table(chry, report)
tab
fisher.test(tab)
# 
#     report
# chry     0     1
#    0 12093    32
#    1     7     6
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  tab
# p-value = 1.051e-12
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#    83.89346 1182.03779
# sample estimates:
# odds ratio 
#    318.899 
