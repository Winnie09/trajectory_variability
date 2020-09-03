setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/testtime/data/data/')
r.tradeseq <- readRDS('tradeSeq/testres.rds')
fdr.tradeseq1 = r.tradeseq[[1]][[1]][,3]
fdr.tradeseq2 = r.tradeseq[[2]][[1]][,3]
fdr.tradeseq3 = r.tradeseq[[3]][[1]][,3]

sum(fdr.tradeseq1 <0.05)
# [1] 3592
sum(fdr.tradeseq2 <0.05)
# [1] 4811
sum(fdr.tradeseq3 <0.05, na.rm = T)
# [1] 1346


r.em = readRDS('EM_SelectKnots/testres.rds')
r.tscan = readRDS('tscan/testres.rds')
r.monocle2 = readRDS('monocle2/testres.rds')
r.monocle3 = readRDS('monocle3/testres.rds')


sum(r.em$fdr < 0.05)
sum(r.tscan$fdr < 0.05)
sum(r.monocle2$fdr < 0.05)
sum(r.monocle3$fdr < 0.05)


perll distribution
gene expression along pseudotime


g = names(sort(r.em$fdr)[1])
perll = r.em$perll
fc = r.em$foldchange
orill = fc + rowMeans(perll)
pt = r.em$pseudotime



setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/')
source('./function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
expr <- readRDS('./hca/data/HCA/proc/matrix/saver.rds')
ct <- readRDS('./hca/data/HCA/proc/ct/sc.rds')
id <- which(ct %in% c('HSC','MEP','Ery'))
expr <- expr[,id]
expr <- expr[rowMeans(expr > 0.1) > 0.1,] ## [1:9070, 1:13269] 
expr = expr[, ]

par(mfrow=c(2,2))
plot(expr[g, names(pt)] ~ pt, pch = 20, xlab = 'pseudotime', ylab = 'expr')
hist(perll[g, ], col = 'grey', breaks = 100)
abline(v = orill[g], col = 'red')


str(expr)
 # num [1:9070, 1:13269] 0.4383 0.0356 0.2666 0.1801 0.1686 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:9070] "DPM1:ENSG00000000419" "C1orf112:ENSG00000000460" "CFH:ENSG00000000971" "FUCA2:ENSG00000001036" ...
 #  ..$ : chr [1:13269] "BM1:52:female_350" "BM1:52:female_463" "BM1:52:female_563" "BM1:52:female_940" ...
head(pt)
# BM4:29:male_144719 BM4:29:male_132013 BM4:29:male_135118  BM4:29:male_60575 
#                  1                  2                  3                  4 
# BM4:29:male_177332 BM4:29:male_109428 
#                  5                  6 
identical(names(pt),colnames(expr))
# [1] FALSE
p <- sub(':.*','',names(pt))
pd=data.frame(e=expr[g,names(pt)],pt=pt,p=p)
head(pd)
                          e pt   p
BM4:29:male_144719 3.077243  1 BM4
BM4:29:male_132013 2.612352  2 BM4
BM4:29:male_135118 2.417650  3 BM4
BM4:29:male_60575  3.305679  4 BM4
BM4:29:male_177332 2.680549  5 BM4
BM4:29:male_109428 3.340420  6 BM4
library(ggplot2)
ggplot(pd,aes(x=pt,y=e)) + geom_point() + facet_wrap(~p)
