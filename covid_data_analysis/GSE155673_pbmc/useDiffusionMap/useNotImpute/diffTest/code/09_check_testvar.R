library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
source(here('function/01_function.R'))

type = 'temra_testvar'
Res.temra <- readRDS(here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type, 'res.rds'))

type = 'tex_testvar'
Res.tex <- readRDS(here('covid','GSE155673_pbmc','useDiffusionMap','useNotImpute','diffTest','result', type, 'res.rds'))


Res.tex$fdr[Res.tex$fdr < 0.05]
# tex_testvar
#   IGKC:ENSG00000211592 JCHAIN:ENSG00000132465  IGLC2:ENSG00000211677 
#           6.707400e-24           1.891695e-02           1.461628e-17 
#  IGLC3:ENSG00000211679 
#           2.001122e-08 
# 
# temra_testvar
#   IGKC:ENSG00000211592 JCHAIN:ENSG00000132465  IGLC2:ENSG00000211677 
#           2.529666e-45           9.233203e-05           1.604591e-08 
#  IGLC3:ENSG00000211679 
#           4.947027e-02 


IGKC:ENSG00000211592 
JCHAIN:ENSG00000132465  
IGLC2:ENSG00000211677 
IGLC3:ENSG00000211679 


fdr.temra <- Res.temra$fdr
fdr.tex <- Res.tex$fdr

head(sort(fdr.temra))
head(sort(fdr.tex))

plotGene(Res.temra, names(sort(fdr.temra)[1:9]))
plotGene(Res.tex, names(sort(fdr.tex)[1:9]))
