library(here)
setwd(here())
source('function/01_function.R')
ddir <- 'harcohen/data/'
rdir <- 'harcohen/testvar/result/'
order <- readRDS(paste0(ddir, 'order.rds'))
meta <- readRDS(paste0(ddir, 'meta.rds'))
expr <- readRDS(paste0(ddir, 'norm.rds'))

m = expr[, order[[2]]]
cellanno <- data.frame(cell = colnames(m),
                       sample = sub(':.*', '', colnames(m)),
                       stringsAsFactors = FALSE)

design = data.frame(intercept = 1,
                    responder = ifelse(meta[,2] == 'Responder', 1, 0),
                    stringsAsFactors = FALSE)
rownames(design) <- meta[,1]

tab = table(cellanno[,2])
select.p <- names(tab[tab > 50])
select.p <- intersect(select.p, meta[meta[,3] == 'anti-PD1', 1])

cellanno <- cellanno[cellanno[,2] %in% select.p, ]
design <- design[select.p, ]
m <- m[, cellanno[,1]]
pseudotime = seq(1, ncol(m))
names(pseudotime) <- colnames(m)
m <- m[rowMeans(m>0.1)>0.01, ]
design <- design[unique(cellanno[,2]), ]
system.time({
  res <- testpt(expr=m, cellanno=cellanno, pseudotime=pseudotime, design=design, type='Variable', ncores = 2, demean = FALSE)
})
saveRDS(res, paste0(rdir, 'branch2_antiPD1_responder_res.rds'))


