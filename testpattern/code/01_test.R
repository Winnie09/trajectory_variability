# ------------
# prepare data
# ------------
library(parallel)
library(splines)
data <- as.character(commandArgs(trailingOnly = TRUE)[[1]][1]) ## clusterType9_1
print(paste0('Analyzing ',data, '...'))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/result/', data, '/')
plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testpattern/plot/', data, '/')
dir.create(rdir, recursive = TRUE)
dir.create(plotdir, recursive = TRUE)
d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/data/data/saver/', data, '.rds'))
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
  Res <- ptest(expr = m, cellanno = ca, pseudotime = pseudotime, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type='Variable', fit.resolution = 1000)
  saveRDS(Res, paste0(rdir, 'ptest_res.rds'))
})

#     user   system  elapsed 
# 74525.76 11940.81  3409.03 

names(Res)
str(Res$res)

# check <<<<<<<<<<<<<
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/data/data/selgene/selgene.rds')
selgene <- intersect(sub(':.*', '', selgene), rownames(res))
apply(res[selgene, ], 2, summary)
apply(res[!rownames(res) %in% selgene, ], 2, summary)

# >>>>>>>>>>>>>>>>

# ----------------------
# plot significant genes
# ----------------------
res <- Res$res
write.csv(Res$res, file = paste0(rdir, 'ptest_res.csv'), quote = FALSE)
# mean diff
res1 <- res[res$meandiff.fdr < 0.05, ]


if (nrow(res1) > 0){
  write.csv(res1, file = paste0(rdir, 'ptest_interceptdiff.csv'), quote = FALSE)
  pdf(paste0(plotdir, 'interceptdiff_fdr_gene.pdf'), width = 10, height = 7)
  gene = rownames(res1)[order(res1[,'meandiff.fdr'])]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res1))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
  
  pdf(paste0(plotdir, 'interceptdiff_diff_gene.pdf'), width = 10, height = 7)
  gene = rownames(res1)[order(abs(res1[,'meandiff.diff']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res1))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
  
  pdf(paste0(plotdir, 'interceptdiff_lfc_gene.pdf'), width = 10, height = 7)
  gene = rownames(res1)[order(abs(res1[,'meandiff.lfc']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res1))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
} else {
  print('No interceptdiff fdr < 0.05!')
}
  
# trend diff
res2 <- res[res$trenddiff.fdr < 0.05, ]
if (nrow(res2) > 0){
  write.csv(res2, file = paste0(rdir, 'ptest_trenddiff.csv'), quote = FALSE)
  pdf(paste0(plotdir, 'trenddiff_fdr_gene.pdf'), width = 10, height = 7)
  gene = rownames(res2)[order(res2[,'trenddiff.fdr'])]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res2))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
  
  pdf(paste0(plotdir, 'trenddiff_diff_gene.pdf'), width = 10, height = 7)
  gene = rownames(res2)[order(abs(res2[,'trenddiff.diff']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res2))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
  
  pdf(paste0(plotdir, 'trenddiff_lfc_gene.pdf'), width = 10, height = 7)
  gene = rownames(res2)[order(abs(res2[,'trenddiff.lfc']), decreasing = TRUE)]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res2))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
} else {
  print('No trenddiff fdr < 0.05!')
}
 
## intercept diff but no trend diff
res3 <- res[res$meandiff.fdr < 0.05 & res$trenddiff.fdr > 0.05, ]
if (nrow(res3) > 0){
  write.csv(res3, file = paste0(rdir, 'ptest_interceptdiff_butNoTrenddiff.csv'), quote = FALSE)
  pdf(paste0(plotdir, 'interceptdiff_butNoTrenddiff_fdr_gene.pdf'), width = 10, height = 7)
  gene = rownames(res3)[order(res3[,'trenddiff.fdr'])]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res3))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
}
  
# ------------------------
# plot insignificant genes
# ------------------------
# intercept diff
res4 <- res[res$meandiff.fdr > 0.05, , drop = FALSE]
if (nrow(res4) > 0){
  write.csv(res4, file = paste0(rdir, 'ptest_interceptdiff_insig.csv'), quote = FALSE)
  pdf(paste0(plotdir, 'interceptdiff_fdr_gene_insig.pdf'), width = 10, height = 7)
  gene = rev(rownames(res4)[order(res4[,'meandiff.fdr'])])
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res4))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
} else {
  print('No interceptdiff fdr > 0.05!')
}
  
# trend diff
res5 <- res[res$trenddiff.fdr > 0.05, , drop = FALSE]
if (nrow(res5) > 0){
  write.csv(res5, file = paste0(rdir, 'ptest_trenddiff_insig.csv'), quote = FALSE)
  pdf(paste0(plotdir, 'trenddiff_fdr_gene_insig.pdf'), width = 10, height = 7)
  gene = rownames(res5)[order(res5[,'trenddiff.fdr'])]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res5))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
} else {
  print('No trenddiff fdr > 0.05!')
}

