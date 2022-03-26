time = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/time.rds')
mem = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/memory.rds')

time[is.na(time[,3]), 3] = 73 ## max time is 3-day limit, use 73 for TIME OUT
mem[is.na(mem[,3]), 3] = 401 ## max mem is 400GB limit, use 410 for OUT OF MEMORY


t2 = tapply(time[,3], list(time[,4], time[,5]), mean, na.rm = T)
m2 = tapply(mem[,3], list(mem[,4], mem[,5]), mean, na.rm = T)
write.csv(t2, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/time_matrix.csv')
write.csv(m2, '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/memory_matrix.csv')

library(reshape2)
t2[is.na(t2)] <- 73
t2 <- round(t2,3)
t2 = melt(t2)
t2[,2] <- factor(as.character(t2[,2]),levels=c('hcaSimu','hca','covid','tb'))


m2[is.na(m2)] <- 401
m2 <- round(m2,3)
m2 <- melt(m2)
m2[,2] <- factor(as.character(m2[,2]),levels=c('hcaSimu','hca','covid','tb'))

library(ggplot2)
library(RColorBrewer)

coldf = read.csv('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/plot/perf/color_code.csv', row.names = 1, as.is = T)
colv = coldf[,2]
names(colv) = coldf[,3]

source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/ggplot_theme.R')
theme_set(.new_theme)

p1 <- ggplot(t2, aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_line(linetype = "dashed", size = 0.2) + 
  geom_point(size = 1) + 
  ylab('time (hour)') + xlab('data') + 
  scale_color_manual(values = colv) 

p2 <- ggplot(m2, aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_line(linetype = "dashed", size = 0.2) + 
  geom_point(size = 1) + 
  ylab('memory (GB)') + xlab('data') + 
  scale_color_manual(values = colv) 

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/efficiency/plot/eff.pdf', width = 6, height = 1.5)
gridExtra::grid.arrange(p1,p2, nrow = 1)
dev.off()



