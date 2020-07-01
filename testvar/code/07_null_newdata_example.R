method = as.character(commandArgs(trailingOnly = T)[[1]])
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testvar/result/newdata/null/' ##########
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
source('./function/01_function.R')
suppressMessages(library(parallel))
suppressMessages(library(splines))
suppressMessages(library(limma))

dir.create(paste0(rdir, method), recursive = T, showWarnings = F)
expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

pseudotime <- readRDS('./testtime/data/data/null/pseudotime.rds')
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)

r.em <- readRDS(paste0(rdir,'EM_SelectKnots/testres.rds'))
g = names(sort(r.em$fdr)[1])
# g = "FUCA2:ENSG00000001036"

r.tradeseq <- readRDS(paste0(rdir,'tradeSeq/testres.rds'))
sum(r.tradeseq[[1]][[1]][,3] < 0.05)
sum(r.tradeseq[[2]][[1]][,3] < 0.05)
sum(r.tradeseq[[3]][[1]][,3] < 0.05, na.rm = TRUE)

# > str(r.tradeseq)
# List of 3
#  $ diffEndTest:List of 1
#   ..$ res:'data.frame': 9070 obs. of  3 variables:
#   .. ..$ waldStat : num [1:9070] 1004 823 606 570 517 ...
#   .. ..$ P.Value  : num [1:9070] 0 0 0 0 0 0 0 0 0 0 ...
#   .. ..$ adj.P.Val: num [1:9070] 0 0 0 0 0 0 0 0 0 0 ...
#  $ patternTest:List of 1
#   ..$ res:'data.frame': 9070 obs. of  3 variables:
#   .. ..$ waldStat : num [1:9070] 16805 16058 9713 7599 7113 ...
#   .. ..$ P.Value  : num [1:9070] 0 0 0 0 0 0 0 0 0 0 ...
#   .. ..$ adj.P.Val: num [1:9070] 0 0 0 0 0 0 0 0 0 0 ...
#  $ earlyDETest:List of 1
#   ..$ res:'data.frame': 9070 obs. of  3 variables:
#   .. ..$ waldStat : num [1:9070] 4842 3502 2945 2652 2550 ...
#   .. ..$ P.Value  : num [1:9070] 0 0 0 0 0 0 0 0 0 0 ...
#   .. ..$ adj.P.Val: num [1:9070] 0 0 0 0 0 0 0 0 0 0 ...
# > sum(r.tradeseq[[1]][[1]][,3] < 0.05)
# [1] 2675
# > sum(r.tradeseq[[2]][[1]][,3] < 0.05)
# [1] 4662
# > sum(r.tradeseq[[3]][[1]][,3] < 0.05, na.rm = TRUE)
# [1] 3077

r.tscan <- readRDS(paste0(rdir,'tscan/testres.rds'))
sum(r.tscan[,3] < 0.05)

r.monocle2 <- readRDS(paste0(rdir, 'monocle2/testres.rds'))
sum(r.monocle2[,3] < 0.05)

r.monocle3 <- readRDS(paste0(rdir,'monocle3/testres.rds'))
sum(r.monocle3[,3] < 0.05)

# > sum(r.tscan[,3] < 0.05)
# [1] 3891
# > sum(r.monocle2[,3] < 0.05)
# [1] 9065
# > sum(r.monocle3[,3] < 0.05)
# [1] 8599

deg.tradeseq1 = rownames(r.tradeseq[[1]][[1]][r.tradeseq[[1]][[1]][, 3]<0.05,])
deg.tradeseq2 = rownames(r.tradeseq[[2]][[1]][r.tradeseq[[1]][[1]][, 3]<0.05,])
deg.tradeseq3 = rownames(r.tradeseq[[3]][[1]][r.tradeseq[[1]][[1]][, 3]<0.05,])
deg.tscan = rownames(r.tscan[r.tscan[,3] < 0.05,])
deg.monocle2 = rownames(r.monocle2[r.monocle2[,3] < 0.05,])
deg.monocle3 = rownames(r.monocle3[r.monocle3[,3] < 0.05,])

v = names(which(table(c(deg.tradeseq1, deg.tradeseq2, deg.tradeseq3, deg.tscan, deg.monocle2, deg.monocle3)) == 6))
g = "AATF:ENSG00000275700"
r.tradeseq[[1]][[1]][g, 3]
r.tradeseq[[2]][[1]][g, 3]
r.tradeseq[[3]][[1]][g, 3]
r.tscan[g, 3]
r.monocle2[g, 3]
r.monocle3[g, 3]
r.em$fdr[g]
plotGene(testptObj = r.em, Gene = g, Mat = expr, Pseudotime = pseudotime, Cellanno = cellanno, Design = design,  Alpha=0.1, Size=0.2, PlotPoints = TRUE, FreeScale = FALSE, BySample = FALSE, type = 'Variable')

