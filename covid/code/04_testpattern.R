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

# -----
# test
# -----
system.time({
  Res <- ptest(expr = m, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(), type='Variable', fit.resolution = 1000)
  saveRDS(Res, paste0(rdir, 'ptest_res.rds'))
})
  
names(Res)
str(Res$res)

# ----------------------
# plot significant genes
# ----------------------
res <- Res$res
# mean diff
res1 <- res[res$meandiff.fdr < 0.05, ]
if (nrow(res1) > 0){
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
  print('No mean diff !')
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
  print('No trend diff!')
}
  
