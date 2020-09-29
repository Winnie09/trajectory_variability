source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/testpt_res.rds')
m <- readRDS('/home-4/zji4@jhu.edu/scratch/diffpt/covid/PRJCA002413_pbmc/data/proc/pt/expr.rds')
pt <- data.frame(Cell = colnames(m), Pseudotime = 1:ncol(m), stringsAsFactors = FALSE)
cellanno = data.frame(cell=colnames(m), sample = sub(':.*','',colnames(m)), stringsAsFactors = FALSE)
rownames(cellanno) <- cellanno[,1]
unis <- unique(cellanno[,2])
design <- cbind(1,as.numeric(grepl('Healthy',unis)))
rownames(design) <- unis
colnames(design) = c(paste0('var_', seq(1, ncol(design))))

## change the gene ID information as only using gene names
rownames(m) <- sub(':.*', '', rownames(m))
names(res$fdr) <- sub(':.*', '', names(res$fdr))
names(res$foldchange) <- sub(':.*', '', names(res$foldchange))
names(res$knotnum) <- sub(':.*', '', names(res$knotnum))
names(res$parameter) <- sub(':.*', '', names(res$parameter))

# -------
# plots
# -------
## heatmap for significant genes
pd <- get_heatmap_plots(testptObj = res, Mat = m, Pseudotime = pt, Cellanno = cellanno, Design = design[,2,drop=FALSE], Max.gene = NULL)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/association_time_heatmap.pdf')
plotQuickHeatmap(pd)
dev.off()

## fitting curves for significant genes
gene <- names(sort(res$fdr)[1:100])
g <- plotGene(testptObj = res, Gene = gene, Mat = m, Pseudotime = pt, Cellanno = cellanno, Design =design[,2,drop = F],  Alpha=1, Size=0.5, PlotPoints = FALSE, FreeScale = T, BySample = FALSE, type ='Variable', colorBySample = FALSE)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/significant_genes_100.pdf', g, dpi = 200, width = 18, height = 18)

## TCF7 as the genes to construct pseudotime (low to high)
g <- plotGene(testptObj = res, Gene = 'TCF7', Mat = m, Pseudotime = pt, Cellanno = cellanno, Design =design[,2,drop = F],  Alpha=1, Size=0.5, PlotPoints = FALSE, FreeScale = T, BySample = FALSE, type ='Variable', colorBySample = FALSE)
ggsave('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/plot/TCF7.pdf', g, dpi = 200, width = 5, height = 5)

# --------------------------------
# save files for significant genes
# --------------------------------
stat = data.frame(Gene_Name = names(res$fdr), FDR = res$fdr, Foldchange = res$foldchange, stringsAsFactors = FALSE)
stat = stat[stat$FDR < 0.05, , drop = FALSE]
stat = stat[order(stat$FDR), , drop = FALSE]
write.csv(stat, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/significant_genes_healthy_vs_covid.csv', row.names = FALSE)
