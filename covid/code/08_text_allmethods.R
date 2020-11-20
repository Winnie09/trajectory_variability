method <- as.character(commandArgs(trailingOnly = T)[[1]])
print(method)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/data/tex/')
cellanno <- readRDS('cellanno.rds')
expr <- readRDS('log2norm.rds')
pt <- readRDS('pseudotime.rds')
ds <- readRDS('design.rds')
design = data.frame(intercept = 1, type = ifelse(ds$type=='HD', 0, 1), stringsAsFactors = F)
rownames(design) <- rownames(ds)
pt <- pt[names(pt) %in% colnames(expr)]
pseudotime <- data.frame(cell = names(pt), time = pt, stringsAsFactors = FALSE)
expr <- expr[rowMeans(expr>0.1)>0.01, ]

source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/result/'
dir.create(rdir, recursive = T)
# testtime
# res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Time', ncores = 8)
# saveRDS(res, paste0(rdir, 'tex_testtime.rds'))
# trenddiff
if (method == 'trenddiff') {
  res <- testpt(expr=expr, cellanno=cellanno, pseudotime=pt, design=design, type='Variable', ncores = 8)
  saveRDS(res, paste0(rdir, 'tex_trenddiff.rds'))
}

## meandiff
if (method == 'meandiff'){
  meanres <- meandiff(expr = expr, cellanno = cellanno, design = design, ncores =8)
  saveRDS(meanres, paste0(rdir, 'tex_meandiff.rds'))
  
}

if (method == 'TSCAN'){
  res <- TSCAN_testvar(expr = expr, cellanno = cellanno, design = design, pseudotime = pseudotime)
  saveRDS(res, paste0(rdir, 'tex_TSCAN.rds'))
}

if (method == 'monocle2'){
  res <- monocle2_testvar(expr = expr, cellanno = cellanno, design = design, pseudotime = pseudotime)
  saveRDS(res, paste0(rdir, 'tex_moncle2.rds'))
}

if (method == 'monocle3'){
  res <- monocle3_testvar(expr = expr, cellanno = cellanno, design = design, pseudotime = pseudotime)
  saveRDS(res, paste0(rdir, 'tex_moncle3.rds'))
}

