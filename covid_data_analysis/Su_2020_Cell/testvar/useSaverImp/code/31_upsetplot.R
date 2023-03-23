rm(list=ls())
library(here)
library(ggplot2)
library(RColorBrewer)
# setwd(here())
# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
setwd('/home/whou10/scratch16/whou10/trajectory_variability/')
pd <- readRDS('covid/Su_2020_Cell/testvar/useSaverImp/plot/perf/pd_heatmap.rds')
pd <- pd[, colnames(pd)!='monocle2_trajtest']
v = colnames(pd)
v[v == "tradeseq_diffEndTest"] <- 'tradeSeqDT'
v[v == "tradeseq_patternTest"] <- 'tradeSeqPT'
v[v == "tradeseq_earlyDETest"] <- 'tradeSeqET'
v[v == "monocle2_trajtest.corr"] <- 'monocle2TrajTestCorr'
colnames(pd) <- v
mat = as.data.frame(pd)
str(mat)
library(UpSetR)
pdf('covid/Su_2020_Cell/testvar/useSaverImp/plot/perf/upset_plot.pdf',width=8,height=3.8)
upset(mat, sets = names(mat), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.off()

