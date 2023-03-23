suppressMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(PseudotimeDE))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(irlba))
print(Sys.time())

dataType=as.numeric(commandArgs(trailingOnly = T))
expr <- readRDS(paste0('/hpc/group/jilab/zj/lamian/pseudotimeDE/count/', dataType, '.rds'))
logexpr <- log2(expr + 1)
cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)
mysce <- SingleCellExperiment(list(counts=expr, logcounts = logexpr)) ## !!
rd <- irlba::prcomp_irlba(t(logcounts(mysce)), scale. = FALSE)$x[, 1:2]

reducedDims(mysce) <- SimpleList(PCA = rd)
colData(mysce)$cl <- 1

fit_ori <- slingshot(mysce, reducedDim = 'PCA', clusterLabels = "cl")
ori_tbl <- tibble(cell = colnames(mysce), pseudotime = rescale(colData(fit_ori)$slingPseudotime_1))

library(parallel)
set.seed(123)
## Set the cores for parallelization. Note that mclapply doesnot work on Windows.
ncores <- 20
print(ncores)
options(mc.cores = ncores)
## Number of subsmaples
n = 100
## Ganerate random subsamples
index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(mysce)[2]), size = 0.8*dim(mysce)[2], replace = FALSE)
})
sub_tbl <- mclapply(index, function(x, sce) {
  sce <- sce[, x]
  rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
  reducedDims(sce) <- SimpleList(PCA = rd)
  
  fit <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "cl")
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, ori_tbl, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}, sce = mysce)

system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = rownames(expr),
                                                 ori.tbl = ori_tbl,
                                                 sub.tbl = sub_tbl, ## To save time, use 100 subsamples
                                                 mat = mysce, ## You can also use a matrix or SeuratObj as the input
                                                 model = "nb",
                                                 mc.cores = ncores))

saveRDS(res, paste0('/hpc/group/jilab/zj/lamian/pseudotimeDE/res/', dataType,'.rds'))  

