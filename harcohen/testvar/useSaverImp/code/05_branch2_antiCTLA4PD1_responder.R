library(here)
setwd(here())
source('function/01_function.R')
ddir <- 'harcohen/data/'
rdir <- 'harcohen/testvar/useSaverImp/result/'
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)

order <- readRDS(paste0(ddir, 'order.rds'))
meta <- readRDS(paste0(ddir, 'meta.rds'))
expr <- readRDS(paste0(ddir, 'saverEachSample/log2_saver.rds'))
c.s <- order[[2]][order[[2]] %in% colnames(expr)]
m = expr[, c.s]

cellanno <- data.frame(cell = colnames(m),
                       sample = sub(':.*', '', colnames(m)),
                       stringsAsFactors = FALSE)

design = data.frame(intercept = 1,
                    responder = ifelse(meta[,2] == 'Responder', 1, 0),
                    stringsAsFactors = FALSE)
rownames(design) <- meta[,1]

tab = table(cellanno[,2])
select.p <- names(tab[tab > 50])
select.p <- intersect(select.p, meta[meta[,3] == 'anti-CTLA4+PD1', 1])

cellanno <- cellanno[cellanno[,2] %in% select.p, ]
design <- design[select.p, ]
m <- m[, cellanno[,1]]
pseudotime = seq(1, ncol(m))
names(pseudotime) <- colnames(m)
m <- m[rowMeans(m>0.1)>0.01, ]
system.time({
  res <- testpt(expr=m, cellanno=cellanno, pseudotime=pseudotime, design=design, type='Variable', ncores = 1, demean = FALSE)
})
saveRDS(res, paste0(rdir, 'branch2_antiCTLA4PD1_responder_res.rds'))


