library(ggplot2)
library(RColorBrewer)
colv = c(brewer.pal(9,'Set1')[c(1:2)], brewer.pal(8, 'Dark2')[1:8])

colv = c(colv, '80bf39', '9cd45d', 'bded87')

allmet = c("Lamian",  "Lamian.chisq", "condiments", "limma",  "monocle2_trajTest",  "monocle2_trajTest.corr", "phenopath", 
 "tradeSeq_diffEndTest", "tradeSeq_earlyDETest", "tradeSeq_patternTest", 'phenopath500', 'phenopath100', 'phenopath300'  )
names(colv) <- c('Lamian', 'Lamian.chisq', setdiff(sort(allmet), c('Lamian', 'Lamian.chisq')))
write.csv(data.frame(method = names(colv), color = colv), '/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/color_code.csv')

