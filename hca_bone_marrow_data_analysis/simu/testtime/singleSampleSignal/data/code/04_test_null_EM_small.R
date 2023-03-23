method <- 'EM_SelectKnnots'
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
rdir <- './testtime/data/data/'
dir.create(paste0(rdir,method), recursive = T, showWarnings = F)
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
allp = sub(':.*','', colnames(expr))
sample <- sub(':.*','',colnames(expr))
names(sample) <- colnames(expr)
design <- cbind(c(1,1,0,0,1,1,0,0))
row.names(design) <- paste0('BM',1:8)
colnames(design) <- 'group'

### reorder cells for each patient
pt <- readRDS(paste0(rdir,'null/pseudotime.rds'))
pt <- data.frame(pt, sample = sub(':.*','', pt$cell), stringsAsFactors = F)
ptlist <- lapply(unique(pt[,3]), function(p){
  tmp <- pt[as.character(pt[,3])==p, ]
  set.seed(12345)
  tmp[, 'pseudotime'] <- sample(tmp[, 'pseudotime'])
  tmp
})
pseudotime <- do.call(rbind, ptlist)
design = cbind(1,design)
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)

# > str(pseudotime)
# 'data.frame': 13269 obs. of  3 variables:
#  $ cell      : chr  "BM4:29:male_164021" "BM4:29:male_280662" "BM4:29:male_76213" "BM4:29:male_56540" ...
#  $ pseudotime: int  250 4362 4437 7898 4975 3665 11342 464 6586 10414 ...
#  $ sample    : chr  "BM4" "BM4" "BM4" "BM4" ...

### 20200701, subset cells and genes
set.seed(1)
id <- sample(colnames(expr),2000)
testres <- testpt(expr=expr[,id],cellanno=cellanno[match(id, cellanno[,1]),],pseudotime=pseudotime[match(id, pseudotime[,1]),],design=design,ncores=4, permuiter=100, type = 'Time')
saveRDS(testres, paste0(rdir, method, '/testres_small.rds'))

