library(Seurat)
library(here)
setwd(here())
data <- readRDS('covid/Su_2020_Cell/data/saver_log2norm.rds')
order <- readRDS('covid/Su_2020_Cell/data/order.rds')
data <- data[, order]
saveRDS(data, 'covid/Su_2020_Cell/data/saver_log2norm_sub.rds')

pt <- seq(1, length(order))
names(pt) <- order
saveRDS(pt, 'covid/Su_2020_Cell/data/tActivate_pseudotime.rds')

meta <- readRDS('covid/Su_2020_Cell/data/meta.rds')
meta <- meta[grepl('Su', meta[,2]), ]
colnames(meta) <- gsub(' ', '_', colnames(meta))
str(meta)
rownames(meta) <- meta$Library_Sample_Code
ap <- sub(':.*','',colnames(data))
meta <- meta[unique(ap), ]
saveRDS(meta, 'covid/Su_2020_Cell/data/meta.rds')

cellanno <- data.frame(cell = colnames(data), sample = ap, stringsAsFactors = FALSE)
saveRDS(cellanno, 'covid/Su_2020_Cell/data/cellanno.rds')

design = data.frame(intercept = 1, type = ifelse(meta$type=='HD', 0, 1), stringsAsFactors = FALSE)
rownames(design) <- rownames(meta)
saveRDS(design, 'covid/Su_2020_Cell/data/design_numeric_HD_others.rds')

design = data.frame(intercept = 1, type = ifelse(meta$type=='HD',0, ifelse(meta$type=='Mi',1, ifelse(meta$type=='Mod', 2, 3))), stringsAsFactors = FALSE)
rownames(design) <- rownames(meta)
saveRDS(design, 'covid/Su_2020_Cell/data/design_numeric_4levels.rds')

design = data.frame(intercept = 1, type = ifelse(meta$type=='HD',0, ifelse(meta$type=='Se',1, NA)), stringsAsFactors = FALSE)
rownames(design) <- rownames(meta)
design <- design[complete.cases(design), ]
saveRDS(design, 'covid/Su_2020_Cell/data/design_numeric_HD_Se.rds')

unitype = unique(meta$type)
for (i in unitype[1:(length(unitype) - 1)]){
  id = which(unitype == i)
  for (j in unitype[(id+1) : length(unitype)]){
    print(paste0(i, ';', j))
    design = data.frame(intercept = 1, type = ifelse(meta$type == i, 0, ifelse(meta$type == j, 1, NA)), stringsAsFactors = FALSE)
    rownames(design) <- rownames(meta)
    design <- design[complete.cases(design), ]
    print(str(design))
    saveRDS(design, paste0('covid/Su_2020_Cell/data/design_numeric_', i, '_', j, '.rds'))
  }
}


a <- meta[meta$type=='Se' & meta[,'Clinical_Outcome(COVID-19)'] %in% c('Recovered','Deceased'),'Library_Sample_Code']
design <- data.frame(intercept = 1, type = ifelse(meta[a, 'Clinical_Outcome(COVID-19)'] == 'Recovered', 0, 1), stringsAsFactors = FALSE)
rownames(design) <- a
saveRDS(design, 'covid/Su_2020_Cell/data/design_numeric_Recovered_Deceased.rds')
