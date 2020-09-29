# ------------
# prepare data
# ------------
library(parallel)
library(splines)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/testvar/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/testvar/'
m <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/PRJCA002413_pbmc/data/proc/pt/expr.rds')
pt <- 1:ncol(m)
names(pt) <- colnames(m)
cellanno = data.frame(cell=colnames(m), sample = sub(':.*','',colnames(m)), stringsAsFactors = FALSE)
rownames(cellanno) <- cellanno[,1]
unis <- unique(cellanno[,2])
design <- cbind(1,as.numeric(grepl('Healthy',unis)))
rownames(design) <- unis
colnames(design) <- c('intercept', 'condition')
m <- m[rowMeans(m > 0.1) > 0.01,]
rownames(m) <- sub(':.*', '', rownames(m))
m <- m[!duplicated(rownames(m)), ]

# # subsample
# m <- m[1:100, ]
# Res <- ptest(expr = m, cellanno = cellanno, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type='Variable', fit.resolution = 1000)

# -----
# test
# -----
system.time({
  Res <- ptest(expr = m, cellanno = cellanno, pseudotime = pt, design=design, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type='Variable', fit.resolution = 1000)
  saveRDS(Res, paste0(rdir, 'ptest_res.rds'))
})
  
names(Res)
sapply(Res, length)
str(Res$res)

# ----------------------
# plot significant genes
# ----------------------
res <- Res$res
write.csv(Res$res, file = paste0(rdir, 'ptest_res.csv'), quote = FALSE)
# intercept diff
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
res4 <- res[res$meandiff.fdr > 0.05, ]
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
res5 <- res[res$trenddiff.fdr > 0.05, ]
if (nrow(res5) > 0){
  write.csv(res5, file = paste0(rdir, 'ptest_trenddiff_insig.csv'), quote = FALSE)
  pdf(paste0(plotdir, 'trenddiff_fdr_gene_insig.pdf'), width = 10, height = 7)
  gene = rownames(res5)[order(res5[,'trenddiff.fdr'])]
  print(plotGene(testptObj = Res, gene = gene[1:min(16, nrow(res5))], variable = 'condition', plot.point = T, point.alpha = 0.1, point.size = 0.1))
  dev.off()
} else {
  print('No trenddiff fdr > 0.05!')
}




