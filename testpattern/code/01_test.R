# ------------
# prepare data
# ------------
library(parallel)
library(splines)
data <- as.character(commandArgs(trailingOnly = TRUE)[[1]][1])
print(paste0('Analyzing ',data, '...'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/result/', data, '/')
plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/plot/', data, '/')
dir.create(rdir, recursive = TRUE)
dir.create(plotdir, recursive = TRUE)
d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/data/data/count/', data, '.rds'))
rownames(d) <- sub(':.*','',rownames(d))
d <- d[!duplicated(rownames(d)), ]
m = log2(d + 1)
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/pseudotime.rds')
pseudotime = pt[,2]
names(pseudotime) = pt[,1]
ap <- sub(':.*', '', colnames(m))
design = cbind(1, c(1,1,0,0,1,1,0,0))
rownames(design) = paste0('BM', seq(1,8))
colnames(design) <- c('intersect', 'condition')
ca <- data.frame(Cell = colnames(m), Sample = ap, stringsAsFactors = FALSE)

# -----
# test
# -----
system.time({
  Res <- ptest(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=4, type='Variable', fit.resolution = 1000)
  saveRDS(Res, paste0(rdir, 'ptest_res.rds'))
})

#     user   system  elapsed 
# 74525.76 11940.81  3409.03 

names(Res)
str(Res$res)

# ----------------------
# plot significant genes
# ----------------------
res <- Res$res
# mean diff
res1 <- res[res$meandiff.fdr < 0.05, ]
if (nrow(res1)){
  pdf(paste0(plotdir, 'meandiff_fdr_gene.pdf'), width = 10, height = 10)
  gene = rownames(res1)[order(res1[,'meandiff.fdr'])]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res1))], variable = 'condition'))
  dev.off()
  
  pdf(paste0(plotdir, 'meandiff_diff_gene.pdf'), width = 10, height = 10)
  gene = rownames(res1)[order(abs(res1[,'meandiff.diff']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res1))], variable = 'condition'))
  dev.off()
  
  pdf(paste0(plotdir, 'meandiff_lfc_gene.pdf'), width = 10, height = 10)
  gene = rownames(res1)[order(abs(res1[,'meandiff.lfc']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res1))], variable = 'condition'))
  dev.off()
} else {
  print('No meandiff fdr < 0.05!')
}
  
# trend diff
res2 <- res[res$trenddiff.fdr < 0.05, ]
if (nrow(res2) > 0){
  pdf(paste0(plotdir, 'trenddiff_fdr_gene.pdf'), width = 10, height = 10)
  gene = rownames(res2)[order(res2[,'trenddiff.fdr'])]
  print(plotGene(testptObj = Res, gene = gene[1:16], variable = 'condition'))
  dev.off()
  
  pdf(paste0(plotdir, 'trenddiff_diff_gene.pdf'), width = 10, height = 10)
  gene = rownames(res2)[order(abs(res2[,'trenddiff.diff']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:16], variable = 'condition'))
  dev.off()
  
  pdf(paste0(plotdir, 'trenddiff_lfc_gene.pdf'), width = 10, height = 10)
  gene = rownames(res2)[order(abs(res2[,'trenddiff.lfc']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:16], variable = 'condition'))
  dev.off()
} else {
  print('No trenddiff fdr < 0.05!')
}
  
