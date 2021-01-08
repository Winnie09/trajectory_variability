library(here)
library(ggplot2)
setwd(here())
source('function/01_function.R')
ddir <- rdir <- 'hca/real/testvar/result/'
pdir <- 'hca/real/testvar/plot/'

### gender
v1 <- v2 <- v1.lm <- v2.lm <- NULL
for (path in c('erythroid', 'lymph', 'monocyte')){
  v1 <- c(v1, readRDS(paste0(rdir, path, '/gender_chrX_overlap.rds')))
  v2 <- c(v2, readRDS(paste0(rdir, path, '/gender_chrY_overlap.rds')))
  v1.lm <- c(v1.lm, readRDS(paste0(rdir, path, '/meandiff_gender_chrX_overlap.rds')))
  v2.lm <- c(v2.lm, readRDS(paste0(rdir, path, '/meandiff_gender_chrY_overlap.rds')))
}

pd <- cbind(ourmethod_chrX = v1, ourmethod_chrY = v2, limma_chrX = v1.lm, limma_chrY = v2.lm)
rownames(pd) <- c('erythroid', 'lymph', 'monocyte')
pd <- reshape2::melt(pd)
pd <- cbind(pd, method = gsub('_.*', '', pd[,2]), chr = gsub('.*_', '', pd[,2]))

library(ggplot2)
pdf(paste0(pdir, path, '/gender_compare_overlap_ourmethod_limma.pdf'), width = 4.5, height = 3)
print(ggplot(data = pd, aes(x = Var2, y = Var1, fill = value)) + geom_tile() + theme_classic() +
  scale_fill_gradient(high = "#132B43",low = "#56B1F7")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab('') + ylab(''))
dev.off()

### age
v1 <-  v1.lm <- NULL
for (path in c('erythroid', 'lymph', 'monocyte')){
  v1 <- c(v1, readRDS(paste0(rdir, path, '/age_overlap.rds')))
  v1.lm <- c(v1.lm, readRDS(paste0(rdir, path, '/meandiff_age_overlap.rds')))
}

pd <- cbind(ourmethod_age = v1, limma_age = v1.lm)
rownames(pd) <- c('erythroid', 'lymph', 'monocyte')
pd <- reshape2::melt(pd)

library(ggplot2)
pdf(paste0(pdir, path, '/age_compare_overlap_ourmethod_limma.pdf'), width = 3, height = 3)
print(ggplot(data = pd, aes(x = Var2, y = Var1, fill = value)) + geom_tile() + theme_classic() +
  scale_fill_gradient(high = "#132B43",
  low = "#56B1F7")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab('') + ylab(''))
dev.off()
