library(ggplot2)
library(RColorBrewer)
colv = c(brewer.pal(9,'Set1')[c(1:2)], brewer.pal(8, 'Dark2')[1:8])

colv = c(colv, '#80BF39', '#9CD45D', '#BDED87')

allmet = c("Lamian",  "Lamian.chisq", "condiments", "limma",  "monocle2_trajTest",  "monocle2_trajTest.corr", "phenopath", 
 "tradeSeq_diffEndTest", "tradeSeq_earlyDETest", "tradeSeq_patternTest", 'phenopath500', 'phenopath100', 'phenopath300')
names(colv) <- c('Lamian', 'Lamian.chisq', setdiff(sort(allmet), c('Lamian', 'Lamian.chisq')))
allmet2 = c('Lamian.pm', 'Lamian.chisq', 'condiments', 'limma', 'monocle2Trajtest', 'monocle2TrajtestCorr', 'phenopath', 'phenopath100', 'phenopath300', 'phenopath500', "tradeSeqDiffEndTest", "tradeSeqEarlyDETest", "tradeSeqPatternTest")
df = data.frame(method = names(colv), color = colv, method2 = allmet2)
write.csv(df, '/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/color_code.csv')


