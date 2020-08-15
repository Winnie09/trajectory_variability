library(parallel)
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

##  test and plot
tmp <- mclapply(c('all', 'start', 'middle', 'end'), function(pos){
      ## slope only
      res <- testpt(expr = m, cellanno = cellanno, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1, type='Variable', test.slope.only = TRUE, test.position = pos)
      saveRDS(res, paste0(rdir, 'slope_', pos, '.rds'))
      gene = names(rev(sort(abs(res$meandiff))))
      pdf(paste0(plotdir, 'slope_', pos, '_meandiff_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      
      gene = names(sort(res$fdr))
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_foldchange.pdf'), width = 4, height = 4)
      print(plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange'))
      dev.off()
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_meandiff.pdf'), width = 4, height = 4)
      print(plot(res$meandiff ~ res$fdr[names(res$meandiff)], pch = 20, xlab = 'fdr', ylab = 'group mean difference'))
      dev.off()
      
      ## all (intersept + slope)
      res <- testpt(expr = m, cellanno = cellanno, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1, type='Variable', test.slope.only = FALSE, test.position = pos)
      saveRDS(res, paste0(rdir, 'all_', pos, '.rds'))
      gene = names(rev(sort(abs(res$meandiff))))
      pdf(paste0(plotdir, 'all_', pos, '_meandiff_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      
      gene = names(sort(res$fdr))
      pdf(paste0(plotdir, 'all_', pos, '_fdr_genes.pdf'), width = 8, height = 8)
      print(plotGene(testptObj = res, gene = gene[1:16], variable = 'condition'))
      dev.off()
      pdf(paste0(plotdir, 'all_', pos, '_fdr_foldchange.pdf'), width = 4, height = 4)
      print(plot(res$foldchange ~ res$fdr[names(res$foldchange)], pch = 20, xlab = 'fdr', ylab = 'LL foldchange'))
      dev.off()
      pdf(paste0(plotdir, 'slope_', pos, '_fdr_meandiff.pdf'), width = 4, height = 4)
      print(plot(res$meandiff ~ res$fdr[names(res$meandiff)], pch = 20, xlab = 'fdr', ylab = 'group mean difference'))
      dev.off()
      return(0)
}, mc.cores = 4)
