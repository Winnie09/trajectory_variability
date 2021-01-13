library(here)
library(reshape2)
library(ggplot2)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

type = 'tex_testtime'
rdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type)
pdir <- here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','plot', type)
# Res <- readRDS(here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type, 'res.rds'))
clu <- readRDS(paste0(rdir, '/cluster.rds'))
fit <- readRDS(paste0(rdir, '/population_fit.rds'))
fit.scale = scalematrix(fit)
fit.scale = fit.scale[names(clu), ]
rownames(fit.scale) <- sub(':.*', '', rownames(fit.scale))

png(paste0(pdir, '/STIP_hm.png'), width = round(ncol(fit.scale)/2), height = round(nrow(fit.scale)*25), res = 200)
print(mySTIP2(fit.scale, rownames(fit.scale)))
dev.off()


