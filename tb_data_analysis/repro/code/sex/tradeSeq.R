library(Matrix)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
method <- 'tradeSeq'

pid <- 2
source('./function/01_function.R')
pt <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/ptpc',pid,'.rds'))
cellanno <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/cellanno.rds')
design <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/design.rds')

cnt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/data/sex/count.rds')

cnt <- cnt[rowMeans(cnt > 0) > 0.01,]

rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/repro/res/sex/pc2/tradeseq/')
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

seed <- as.numeric(commandArgs(trailingOnly = T)[1])
partid <- as.numeric(commandArgs(trailingOnly = T)[2])
set.seed(seed)
samp1 <- sample(rownames(design),nrow(design)/2)
samp2 <- setdiff(rownames(design),samp1)

cellanno1 <- cellanno[cellanno[,2] %in% samp1,]
cellanno2 <- cellanno[cellanno[,2] %in% samp2,]

### two group along pseudotime
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(slingshot))
suppressMessages(library(tradeSeq))
### prepare data [Note: in two group senario, these should be different: design, v, cellWeights]


if (partid==1) {
  
  pdt <- data.frame(curve1 = pt[cellanno1[,1]], curve2 = pt[cellanno1[,1]])
  rownames(pdt) <- cellanno1[,1]
  
  v <- (cellanno1$sample %in% rownames(design)[design[,2]==0] + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- cellanno1[,1]
  ### run test
  set.seed(12345)
  sce <- fitGAM(counts = cnt[,cellanno1[,1]], pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE,BPPARAM=BiocParallel::SerialParam())
  print(str('sce'))
  
  Final <- list()
  for (TestType in c('diffEndTest', 'patternTest', 'earlyDETest')){
    print(TestType)
    if (grepl('diffEndTest', TestType)){
      Res <- diffEndTest(sce)  
    } else if (grepl('patternTest', TestType)){
      Res <- patternTest(sce)  
    } else if (grepl('earlyDETest', TestType)){
      Res <- earlyDETest(sce, knots = c(1,2), global = TRUE, pairwise = TRUE)
    }
    res <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
    row.names(res) <- row.names(Res)
    res <- res[order(res[,3], -res[,1]), ]
    Final[[TestType]] <- res
  }
  
} else {
  pdt <- data.frame(curve1 = pt[cellanno2[,1]], curve2 = pt[cellanno2[,1]])
  rownames(pdt) <- cellanno2[,1]
  
  v <- (cellanno2$sample %in% rownames(design)[design[,2]==0] + 0)
  v <- ifelse(v==1, 0.99, 0.01)
  
  cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
  rownames(cellWeights) <- cellanno2[,1]
  ### run test
  set.seed(12345)
  sce <- fitGAM(counts = cnt[,cellanno2[,1]], pseudotime = pdt, cellWeights = cellWeights,
                nknots = 6, verbose = FALSE,parallel=TRUE,BPPARAM=BiocParallel::SerialParam())
  print('sce')
  
  Final <- list()
  for (TestType in c('diffEndTest', 'patternTest', 'earlyDETest')){
    print(TestType)
    if (grepl('diffEndTest', TestType)){
      Res <- diffEndTest(sce)
    } else if (grepl('patternTest', TestType)){
      Res <- patternTest(sce)
    } else if (grepl('earlyDETest', TestType)){
      Res <- earlyDETest(sce, knots = c(1,2), global = TRUE, pairwise = TRUE)
    }
    res <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
    row.names(res) <- row.names(Res)
    res <- res[order(res[,3], -res[,1]), ]
    Final[[TestType]] <- res
  }
}

saveRDS(Final, paste0(rdir, seed,'_',partid,'.rds'))  


