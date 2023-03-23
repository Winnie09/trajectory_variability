library(reshape2)
library(ggplot2)
####################################################################################3
gsmethod <- 'tradeSeq'
compareid <- (1:100)*10
perf <- list()
method = 'EM_pm'
perf[[method]] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1$statistics
  s1 <- s1[order(s1[,1],-s1[,3]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[1]])[which(s2[[1]]$adj.P.Val < 0.05)]
  # median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
})

method = 'tscan'
perf[[method]] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[order(s1[,3],-s1[,2]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[1]])[which(s2[[1]]$adj.P.Val < 0.05)]
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 

method = 'monocle2'
perf[[method]] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[order(s1[,3],-s1[,2]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[1]])[which(s2[[1]]$adj.P.Val < 0.05)]
  
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 

method = 'tradeSeq'
perf[['tradeSeq_startVsEndTest']] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[['startVsEndTest']]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[1]])[which(s2[[1]]$adj.P.Val < 0.05)]
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
})


method = 'tradeSeq'
perf[['tradeSeq_associationTest']] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[['associationTest']]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[1]])[which(s2[[1]]$adj.P.Val < 0.05)]
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
})

pd <- melt(do.call(rbind,perf))
p1 <- ggplot(pd,aes(x=Var1,y=value)) + geom_boxplot() + theme_classic() + ggtitle('tradeSeq_startend')


####################################################################################3
gsmethod <- 'tradeSeq'
perf <- list()
method = 'EM_pm'
perf[[method]] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1$statistics
  s1 <- s1[order(s1[,1],-s1[,3]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[2]])[which(s2[[2]]$adj.P.Val < 0.05)]
  #median(which(rownames(s1) %in% gs))
  
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 


method = 'tscan'
perf[[method]] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[order(s1[,3],-s1[,2]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[2]])[which(s2[[2]]$adj.P.Val < 0.05)]
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 

method = 'monocle2'
perf[[method]] <- sapply(c(3,5,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[order(s1[,3],-s1[,2]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  #gs <- rownames(s2)[s2$fdr < 0.05]
  gs <- rownames(s2[[2]])[which(s2[[2]]$adj.P.Val < 0.05)]
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 

pd <- melt(do.call(rbind,perf))
p2 <- ggplot(pd,aes(x=Var1,y=value)) + geom_boxplot() + theme_classic() + ggtitle('tradeSeq_assocation')



####################################################################################3
gsmethod <- 'tscan'
perf <- list()
method = 'EM_pm'
perf[[method]] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1$statistics
  s1 <- s1[order(s1[,1],-s1[,3]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  #median(which(rownames(s1) %in% gs))
  
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 


method = 'tradeSeq'
perf[['tradeSeq_startend']] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[[1]]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 

perf[['tradeSeq_association']] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[[2]]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 


method = 'monocle2'
perf[[method]] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[order(s1[,3],-s1[,2]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 

pd <- melt(do.call(rbind,perf))
p3 <- ggplot(pd,aes(x=Var1,y=value)) + geom_boxplot() + theme_classic() + ggtitle('tscan')


####################################################################################3
gsmethod <- 'monocle2'
perf <- list()
method = 'EM_pm'
perf[[method]] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1$statistics
  s1 <- s1[order(s1[,1],-s1[,3]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  #median(which(rownames(s1) %in% gs))
  
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 


method = 'tradeSeq'
perf[['tradeSeq_startend']] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[[1]]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 

perf[['tradeSeq_association']] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[[2]]
  s1 <- s1[order(s1[,3],-s1[,1]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 


method = 'tscan'
perf[[method]] <- sapply(seq(1,8),function(lvid) {
  s1 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',method,'/monocyte/leaveout_',lvid,'/testtime_res_7sample.rds'))
  s1 <- s1[order(s1[,3],-s1[,2]),]
  s2 <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/result/',gsmethod,'/monocyte/leaveout_',lvid,'/testtime_res_1sample.rds'))
  gs <- rownames(s2)[s2$fdr < 0.05]
  #median(which(rownames(s1) %in% gs))
  mean(sapply(compareid, function(i){
    length(intersect(rownames(s1)[1:i], gs))/i
  }))
}) 


pd <- melt(do.call(rbind,perf))
p4 <- ggplot(pd,aes(x=Var1,y=value)) + geom_boxplot() + theme_classic() + ggtitle('monocle2')
library(gridExtra)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime_leave_one_out/plot/overlap_all.pdf', width = 12, height = 9)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()



