# setwd('/Users/wenpinhou/Dropbox/trajectory_variability/')
library(here())
setwd(here())
time = readRDS('efficiency/plot/time.rds')
mem = readRDS('efficiency/plot/memory.rds')

time[is.na(time[,3]), 3] = 80 ## max time is 3-day limit, use 80 for TIME OUT
mem[is.na(mem[,3]), 3] = 450 ## max mem is 400GB limit, use 450 for OUT OF MEMORY


t2 = tapply(time[,3], list(time[,4], time[,5]), mean, na.rm = T)
m2 = tapply(mem[,3], list(mem[,4], mem[,5]), mean, na.rm = T)
colnames(t2) <- colnames(m2) <- c('COVID', 'HCA', 'HCA.Simu', 'TB')
write.csv(t2, 'efficiency/plot/time_matrix.csv')
write.csv(m2, 'efficiency/plot/memory_matrix.csv')

write.csv(t(t2), 'efficiency/plot/time_matrix_transpose.csv')
write.csv(t(m2), 'efficiency/plot/memory_matrix_transpose.csv')


library(reshape2)
t2[is.na(t2)] <- 80
t2 <- round(t2,3)
t2 = melt(t2)
t2[,2] <- factor(as.character(t2[,2]),levels=c('HCA.Simu','HCA','COVID','TB'))


m2[is.na(m2)] <- 450
m2 <- round(m2,3)
m2 <- melt(m2)
m2[,2] <- factor(as.character(m2[,2]),levels=c('HCA.Simu','HCA','COVID','TB'))

library(ggplot2)
library(RColorBrewer)

coldf = read.csv('hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/color_code.csv', row.names = 1, as.is = T)
colv = coldf[,2]
names(colv) = coldf[,3]
names(colv)[names(colv) == 'tradeSeqPatternTest'] = 'tradeSeq'

# source('/Users/wenpinhou/Dropbox/resource/ggplot_theme.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)

p1 <- ggplot(t2, aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_line(linetype = "dashed", size = 0.2) + 
  geom_point(size = 1.2, alpha = 0.8, stroke = 0) + 
  geom_hline(yintercept = 72, linetype = 'dotdash', color = 'black', size = 0.4) + 
  ylab('time (hour)') + xlab('dataset') + 
  scale_color_manual(values = colv) +
  scale_y_continuous(breaks=c(0, 20, 40, 60, 72, 80),
        labels=c("0", '20', '40', '60', '72', 'timeout'))

p2 <- ggplot(m2, aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_line(linetype = "dashed", size = 0.2) + 
  geom_point(size = 1.2, alpha = 0.8, stroke = 0) + 
  geom_hline(yintercept = 400, linetype = 'dotdash', color = 'black', size = 0.4) + 
  ylab('memory (GB)') + xlab('dataset') + 
  scale_color_manual(values = colv) +
  scale_y_continuous(breaks=c(seq(0,4)*100, 440, 465),
        labels=c("0", '100', '200', '300', '400', 'memory', 'out of')) 

pdf('efficiency/plot/eff.pdf', width = 6.5, height = 1.5)
gridExtra::grid.arrange(p1,p2, nrow = 1)
dev.off()



