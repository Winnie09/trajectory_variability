library(parallel)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testvar/data/data/count/clusterType9_4.rds')
rownames(d) <- sub(':.*','',rownames(d))
m <- d[1:100, ]
m = log2(m + 1)
pt <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/null/pseudotime.rds')
ap <- sub(':.*', '', colnames(m))
design = cbind(1, c(1,1,0,0,1,1,0,0))
rownames(design) = paste0('BM', seq(1,8))
colnames(design) <- c('intersect', 'condition')
ca <- data.frame(Cell = colnames(m), Sample = ap)

# expr = m
# cellanno = ca
# pseudotime = pt 
# permuiter=100
# EMmaxiter=100
# EMitercutoff=1
# verbose=F
# ncores=detectCores()
# type='Time'

##  1
res <- testpt(expr = m[1:10,], cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'slope', test.position = 'all')
slope_allpos = names(sort(res$fdr))
{
# plotHeatmap(testptObj = res, Mat = m, Pseudotime = pt, Cellanno = ca, Design = design, Max.gene = 500, heatmap.x = 0.4, heatmap.y = 0.5, heatmap.width = 0.8, heatmap.height = 1.0, dend.x = 0.90, dend.y = 0.44, dend.width = 0.2, dend.height = 0.98)
# pd <- get_heatmap_plots(testptObj = res, Mat = m, Pseudotime = pt, Cellanno = ca, Design = design[,2,drop=F], Max.gene = NULL)
# plotQuickHeatmap(pd)
##  2
res <- testpt(expr = m, cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'slope', test.position = 'start')
slope_start = names(sort(res$fdr))

## 3
res <- testpt(expr = m, cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'slope', test.position = 'middle')
slope_middle = names(sort(res$fdr)[1:9])
## 4
res <- testpt(expr = m, cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'slope', test.position = 'end')
slope_end = names(sort(res$fdr)[1:9])
##  21
res <- testpt(expr = m, cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'all', test.position = 'all')
all_allpos = names(sort(res$fdr))
# plotHeatmap(testptObj = res, Mat = m, Pseudotime = pt, Cellanno = ca, Design = design, Max.gene = 500, heatmap.x = 0.4, heatmap.y = 0.5, heatmap.width = 0.8, heatmap.height = 1.0, dend.x = 0.90, dend.y = 0.44, dend.width = 0.2, dend.height = 0.98)
# pd <- get_heatmap_plots(testptObj = res, Mat = m, Pseudotime = pt, Cellanno = ca, Design = design[,2,drop=F], Max.gene = NULL)
# plotQuickHeatmap(pd)
##  22
res <- testpt(expr = m, cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'all', test.position = 'start')
all_start = names(sort(res$fdr)[1:9])

## 23
res <- testpt(expr = m, cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'all', test.position = 'middle')
all_middle = names(sort(res$fdr)[1:9])
## 24
res <- testpt(expr = m, cellanno = ca, pseudotime = pt, design=design, permuiter=10, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=detectCores(),type='Variable', test.type = 'all', test.position = 'end')
all_end = names(sort(res$fdr))

gene = c(slope_allpos[1:7], slope_start[1:7], slope_middle[1:7], slope_end[1:7], all_allpos[1:7], all_start[1:7], all_middle[1:7], all_end[1:7])

plotGene(testptObj = res, Gene = slope_allpos[1:9], Mat = m, Pseudotime = pt, Cellanno = ca, Design = design[,2, drop = F],  Alpha=1, Size=0.5, PlotPoints = T, FreeScale = TRUE, BySample = FALSE, type = 'Variable', colorBySample = FALSE)
}
plotGene(testptObj = res, Gene = slope_allpos[1], Mat = m, Pseudotime = pt, Cellanno = ca, Design = design[,2, drop = F],  Alpha=1, Size=0.5, PlotPoints = T, FreeScale = TRUE, BySample = FALSE, type = 'Variable', colorBySample = FALSE)


testptObj = res
Gene = slope_allpos[1:3]
Mat = m
Pseudotime = pt
Cellanno = ca
Design = design[,2,drop = FALSE]
Alpha=1
Size=0.5
PlotPoints = T
FreeScale = FALSE
BySample = FALSE
type = 'Variable'
colorBySample = FALSE

# ---------------------------
# check these files: 20200724
# ---------------------------
fitpt_beta.R
fitpt.R
plotGene_Variable_coloredByVariable_v3.R
predict_v3.R
testpt_v3.R

