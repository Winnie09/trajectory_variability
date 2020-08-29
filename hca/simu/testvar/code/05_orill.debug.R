res.null <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/result/null/EM_SelectKnots/testres.rds')
res.em <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/result/EM_SelectKnots/clusterType2_4_testres.rds')

g = c("CENPBD1P1:ENSG00000213753",  "HOMEZ:ENSG00000215271", "RPL37P6:ENSG00000241431", "MESP1:ENSG00000166823", "BRF1:ENSG00000185024", "AC020915.1:ENSG00000142396")

perll.null = res.null[[2]][g, ]
perll.em <- res.em[[2]][g,]

# par(mfrow = c(2,6))
# for (i in g){
#   hist(perll.null[i, ], main = paste0(i, ' null'), breaks = 100, col = 'grey')
#   orill = res.null$foldchange[i] + mean(perll.null[i, ])
#   abline(v = orill, col = 'red')
# }
# 
# for (i in g){
#   hist(perll.em[i, ], main = paste0(i, ' em'), breaks = 100, col = 'grey')
#   orill = res.em$foldchange[i] + mean(perll.em[i, ])
#   abline(v = orill, col = 'red')
# }

# want to know why the same genes are not added signals, but in null it is not differential, and in added signal matrix it is called differential. 
# in null, and added matrix seperately,find the orill
# read added matrix
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('./function/01_function.R')
selgene <- readRDS('./testvar/data/data/selgene/selgene.rds')
ddir <- './testvar/data/data/'
expr <- readRDS(paste0(ddir, 'saver/clusterType', 2, '_', 4, '.rds'))
design = matrix(c(1,1,0,0,1,1,0,0), nrow=8)
dimnames(design) = list(paste0('BM',seq(1,8)), 'group')
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
expr <- log2(expr + 1)
design = cbind(1,design)
pseudotime <- readRDS('./testtime/data/data/null/pseudotime.rds')
pseudotime <- data.frame(Cell = pseudotime[,1], Pseudotime = as.numeric(pseudotime[,2]), stringsAsFactors = FALSE)
pseudotime <- pseudotime[order(pseudotime[,2]), ]
psn <- pseudotime[,2]
names(psn) <- pseudotime[,1]
pseudotime <- psn
# run test
set.seed(12345)
orifit <- fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=2,parallel=FALSE)
knotnum <- orifit$knotnum
orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)
orill.em <- orill[row.names(expr)]

# read null matrix
expr.null <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr.null <- expr.null[,id]
expr.null <- expr.null[rowMeans(expr.null > 0.1) > 0.1,] ## [1:9070, 1:13269] 
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'
# run test
set.seed(12345)
orifit <- fitpt(expr=expr.null, cellanno=cellanno, pseudotime=pseudotime, design=design, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=2,parallel=FALSE)
knotnum <- orifit$knotnum
orill <- sapply(orifit$parameter,function(i) unname(i$ll),USE.NAMES = F)
orill.null <- orill[row.names(expr)]
saveRDS(list(orill.em = orill.em, orill.null = orill.null), '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/result/orill_debug.rds')

