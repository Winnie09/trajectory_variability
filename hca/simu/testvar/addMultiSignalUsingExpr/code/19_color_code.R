pdir = '/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/'
library(ggplot2)
library(RColorBrewer)
colv = c(brewer.pal(9,'Set1')[c(1:2)], brewer.pal(8, 'Dark2')[1:8])
allmet = c("Lamian",  "Lamian.chisq", "condiments", "limma",  "monocle2_trajTest",  "monocle2_trajTest.corr", "phenopath", 
 "tradeSeq_diffEndTest", "tradeSeq_earlyDETest", "tradeSeq_patternTest"  )
names(colv) <- c('Lamian', 'Lamian.chisq', setdiff(sort(allmet), c('Lamian', 'Lamian.chisq')))
write.csv(data.frame(method = names(colv), color = colv), paste0(pdir, 'color_code.csv'))

