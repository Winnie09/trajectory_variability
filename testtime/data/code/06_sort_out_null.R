r.tradeseq <- readRDS('testres.rds')
fdr.tradeseq1 = r.tradeSeq[[1]][[1]][,3]
fdr.tradeseq2 = r.tradeSeq[[2]][[1]][,3]
fdr.tradeseq3 = r.tradeSeq[[3]][[1]][,3]

# > sum(fdr.tradeseq1 <0.05)
# [1] 3592
# > sum(fdr.tradeseq2 <0.05)
# [1] 4811
# > sum(fdr.tradeseq3 <0.05, na.rm = T)
# [1] 1346



