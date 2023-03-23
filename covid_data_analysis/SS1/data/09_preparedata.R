library(Seurat)
library(here)
setwd(here())
data <- readRDS('covid/SS1/data/saver_log2norm.rds')
order <- readRDS('covid/SS1/data/order.rds')
data <- data[, order]
saveRDS(data, 'covid/SS1/data/saver_log2norm_sub.rds')

pt <- seq(1, length(order))
names(pt) <- order
saveRDS(pt, 'covid/SS1/data/tActivate_pseudotime.rds')

meta <- readRDS('covid/Su_2020_Cell/data/metaall.rds')
colnames(meta) <- gsub(' ', '_', colnames(meta))
str(meta)
rownames(meta) <- meta$Library_Sample_Code
ap <- sub(':.*','',colnames(data))
meta <- meta[unique(ap), ]

cellanno <- data.frame(cell = colnames(data), sample = ap, stringsAsFactors = FALSE)
saveRDS(cellanno, 'covid/SS1/data/cellanno.rds')

unitype = unique(meta$type)
for (i in unitype[1:(length(unitype) - 1)]){
  id = which(unitype == i)
  for (j in unitype[(id+1) : length(unitype)]){
    print(paste0(i, ';', j))
    design = data.frame(intercept = 1, type = ifelse(meta$type == i, 0, ifelse(meta$type == j, 1, NA)), stringsAsFactors = FALSE)
    rownames(design) <- rownames(meta)
    design <- design[complete.cases(design), ]
    print(str(design))
    saveRDS(design, paste0('covid/SS1/data/design_numeric_', i, '_', j, '.rds'))
  }
}

