# new a function
# input: testpt output including beta, phi
# input: covariables values user wants to know , if NULL then the unique values of the covatiates in the testpt data. if in the data only have age = 10, 20, 30, users can input 25 then we can output the pseudotime pattern of age == 25.
# phi * x * beta

signal = 1
method = 'EM'
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/perf/'
selgene <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/simu/testvar/addMultiSignalUsingExpr/selgene/selgene.rds')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca')
rdir <- './simu/testvar/addMultiSignalUsingExpr/result/'
ddir <- './data/simu/testvar/addMultiSignalUsingExpr/'
r <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/result/EM/1.rds')
res = data.frame(fdr = r$fdr, foldchange = r$foldchange)
expr = r$expression
pseudotime = r$pseudotime
cellanno = r$cellanno

suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))
suppressMessages(library(RColorBrewer))

res <- res[order(res$fdr),]
tg <- rownames(res)[res$fdr < 0.05]
p <- sub(':.*','',colnames(expr))
expr <- expr[,order(p,pseudotime[colnames(expr)])]

set.seed(12345)
clu <- kmeans(expr[tg,],10,iter.max = 1000)$cluster
clu <- sort(clu)
saveRDS(clu, paste0(rdir, 'clu/geneclu.rds'))

colann <- data.frame(sample=sub(':.*','',colnames(expr)))
rownames(colann) <- colnames(expr)
colann$group <- as.factor(ifelse(colann[,1] %in% paste0('BM', c(3:4, 7:9)), 'Addsignal', 'Nosignal'))

rowann <- data.frame(cluster = clu[tg],
                     spikein = ifelse(tg %in% selgene, 'Yes', 'No'))
rownames(rowann) <- tg

library(pheatmap)
c1 <- brewer.pal(8,'Set1')
names(c1) <- paste0('BM',1:8)
c2 <- brewer.pal(8,'Dark2')[1:2]
names(c2) <- unique(colann$group)
r1 <- brewer.pal(10,'Set3')
names(r1) <- unique(rowann$cluster)
cpl <- colorRampPalette(rev(brewer.pal(n = 7, name =  "RdYlBu")))(100)
cpl <- c(rep(cpl[1],30),cpl,rep(cpl[length(cpl)],10))

# gaps_row=cumsum(rle(clu[1:30])$lengths)
# gaps_col = cumsum(rle(as.character(colann[,1]))$lengths)

png('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/clu/geneclu.30.10.png')
# png('geneclu.30.10.png')
pheatmap(expr[names(clu),],
              color=cpl, 
              cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
              annotation_col = colann,
              annotation_row = rowann[names(clu),,drop = F],
              annotation_colors = list(patient=c1, group = c2, cluster = r1))
dev.off()

