rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
n.permute <- 10000 ## update to be 1e4
max.clunum <- 50
source("/Users/wenpinhou/Dropbox/trajectory_variability/function/01_function.R")
plotdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/auto_pc_auto_nclu_module_3traj/plot/'
rdir <- '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/auto_pc_auto_nclu_module_3traj/result/'
# --------------------------------------------------------------
# input: seurat integrated object including:
# low dim reduction: umap, pca, or phate
# celltype: a dataframe, col 1 is cell name, col 2 is cell type (at least for the cells with origin cell type), col 3 is sample name
# origin: the origin cell type
# --------------------------------------------------------------
# read in data
umap = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/integrate/ser/umap.rds')
pca <- as.matrix(umap@reductions$pca@cell.embeddings)
a = readRDS('/Users/wenpinhou/Dropbox/trajectory_variability/hca/data/HCA/proc/ct/sc.rds')
ct = data.frame(cell = names(a), celltype = a, sample = sapply(names(a), function(i) sub(':.*', '', i)), stringsAsFactors = FALSE)

# load HCA-BM tree
res <- readRDS(paste0(rdir, 'infer_tree_structure_res.rds'))

## ======================
## branch proportion test
## ======================
## setting
source('/Users/wenpinhou/Dropbox/trajectory_variability/function/evaluate_uncertainty.R')
source('/Users/wenpinhou/Dropbox/trajectory_variability/package/Lamian/R/branchPropTest.R')  
design = data.frame(intercetp = rep(1,8), reduce = rep(c(1,1,0,0),2))
rownames(design) <- paste0('BM', seq(1,8))


## perform branch proportion test on HCA-BM data (original values)
tmp.tb <- sapply(res[['order']], function(o){
  table(sapply(o, function(i) sub(':.*','',i)))
})
invisible(capture.output(res.bp <- branchPropTest(data = tmp.tb, design, method = 'multinom')))
# > res.bp
#               (Intercept) covariate
# branch: 5,1     0.7396051 0.2021713
# branch: 5,3,4   0.4952636 0.3861338


## ----------------------------
## for some samples: BM1,2,5,6
## subsample cells, and then redo infer tree structure
## calculate branch composition and save the object (ctcomplist)
## ----------------------------
cellprop.d = list()

for (rm.perc in seq(0.2, 0.8, 0.1)){
  print(rm.perc)
  plotdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256_multinomCellComp/', rm.perc, '/plot/')
  rdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256_multinomCellComp/', rm.perc, '/result/')
  dir.create(plotdir, recursive = T, showWarnings = F)
  dir.create(rdir, recursive = T, showWarnings = F)
  
  selectcell = res$order[[2]]
  selectcell = selectcell[ct[selectcell, 'sample'] %in% c('BM1', 'BM2', 'BM5', 'BM6')]
  
  ## construct pseudotime tree in each simulation
  ## and calculate branch proportion
  ctcomplist <- reproduce.js <- reproduce.oc <- corr.score <- list()
  for (pmid in 1:1e2){
    set.seed(pmid)
    rmcell = sample(selectcell, rm.perc*length(selectcell))
    subset.cell = setdiff(rownames(pca), rmcell)
    print(pmid)
    inferobj = res
    pr.pm <- inferobj$pca[subset.cell,]
    newbranch <- inferobj$branch
    js.cut <- inferobj$js.cut
    oc.cut <- inferobj$oc.cut
    pt <- inferobj$pseudotime
    ord <- inferobj$order
    alls <- inferobj$allsample
    
    ## cluster cells
    invisible(capture.output(clu <-
                               mykmeans(pr.pm, number.cluster = max(inferobj$clusterid))$cluster))
    
    
    ## build pseudotime
    mcl.pm <-
      exprmclust(t(pr.pm), cluster = clu, reduce = FALSE) ###
    
    ## select origin cluster
    pt.pm.mean <-
      tapply(pt[names(mcl.pm[['clusterid']])], list(mcl.pm[['clusterid']]), mean)
    start.cluster <- names(which.min(pt.pm.mean))
    
    ## construct pseudotime
    ord.pm <-
      TSCANorder(
        mcl.pm,
        startcluster = start.cluster,
        listbranch = TRUE,
        orderonly = TRUE
      )
    
    pt.pm <-
      unlist(sapply(sapply(ord.pm, length), function(i)
        seq(1, i)))
    names(pt.pm) <- unname(unlist(ord.pm))
    
    
    # get candidate branches
    newbranch.pm <-
      findbranch(mst = mcl.pm$MSTtree,
                 order = ord.pm,
                 origin = start.cluster)
    
    ## compare two MST
    js <- sapply(seq(1, length(newbranch)), function(i) {
      id <-
        which(sapply(paste0(names(ord), ','), function(k)
          grepl(paste0(
            paste0(newbranch[[i]], collapse = ','), ','
          ), k)))[1]
      cells <- ord[[id]]
      b.ori <-
        intersect(unlist(sapply(newbranch[[i]], function(k)
          names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
      sapply(seq(1, length(newbranch.pm)), function(j) {
        id <-
          which(sapply(paste0(names(ord.pm), ','), function(k)
            grepl(paste0(
              paste0(newbranch.pm[[j]], collapse = ','), ','
            ), k)))[1]
        cells <- ord.pm[[id]]
        b.pm <-
          intersect(unlist(sapply(newbranch.pm[[j]], function(k)
            names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
        js <-
          length(intersect(b.pm, b.ori)) / length(union(b.pm, b.ori))
      })
    })
    oc <- sapply(seq(1, length(newbranch)), function(i) {
      id <-
        which(sapply(paste0(names(ord), ','), function(k)
          grepl(paste0(
            paste0(newbranch[[i]], collapse = ','), ','
          ), k)))[1]
      cells <- ord[[id]]
      b.ori <-
        intersect(unlist(sapply(newbranch[[i]], function(k)
          names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
      sapply(seq(1, length(newbranch.pm)), function(j) {
        id <-
          which(sapply(paste0(names(ord.pm), ','), function(k)
            grepl(paste0(
              paste0(newbranch.pm[[j]], collapse = ','), ','
            ), k)))[1]
        cells <- ord.pm[[id]]
        b.pm <-
          intersect(unlist(sapply(newbranch.pm[[j]], function(k)
            names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
        oc <-
          length(intersect(b.pm, b.ori)) / min(length(b.pm), length(b.ori))
      })
    })
    corr <- sapply(seq(1, length(newbranch)), function(i) {
      id <-
        which(sapply(paste0(names(ord), ','), function(k)
          grepl(paste0(
            paste0(newbranch[[i]], collapse = ','), ','
          ), k)))[1]
      cells <- ord[[id]]
      b.ori <-
        intersect(unlist(sapply(newbranch[[i]], function(k)
          names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
      
      sapply(seq(1, length(newbranch.pm)), function(j) {
        id <-
          which(sapply(paste0(names(ord.pm), ','), function(k)
            grepl(paste0(
              paste0(newbranch.pm[[j]], collapse = ','), ','
            ), k)))[1]
        cells <- ord.pm[[id]]
        b.pm <-
          intersect(unlist(sapply(newbranch.pm[[j]], function(k)
            names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
        ov = intersect(b.ori, b.pm)
        cor(pt[ov], pt.pm[ov])
      })
    })
    corr[is.na(corr)] <- 0
    colnames(corr) <-
      colnames(oc) <-
      colnames(js) <- paste0('original', seq(1, length(newbranch)))
    
    ## get js binary to match branches
    js.binary <- get_binary(js, js.cut)
    corr.score[[pmid]] <- corr * js.binary
    js.melt <- melt(js.binary)
    js.melt <- js.melt[js.melt[, 3] != 0, ]
    colnames(js.melt) <-
      c('permutation.branch', 'original.branch', 'matched')
    reproduce.js[[pmid]] <- as.character(js.melt[, 2])
    
    ## get oc binary to match branches
    oc.binary <- get_binary(oc, oc.cut)
    oc.melt <- melt(oc.binary)
    oc.melt <- oc.melt[oc.melt[, 3] != 0, ]
    reproduce.oc[[pmid]] <- as.character(oc.melt[, 2])
    
    ## samples cell compositions
    ctcomp.new.logit <- ctcomp.new <-
      matrix(0, nrow = length(unique(alls)), ncol = length(newbranch))
    colnames(ctcomp.new.logit) <- colnames(ctcomp.new) <-
      paste0('origin', seq(1, length(newbranch)))
    rownames(ctcomp.new.logit) <- rownames(ctcomp.new) <- unique(alls)
    
    if (nrow(js.melt) > 0) {
      ctcomp <-
        sapply(seq_len(nrow(js.melt)), function(i) {
          ## corrected  from 2 to 1. 2020/08/31
          c <- names(clu)[clu %in% newbranch.pm[[js.melt[i, 1]]]]
          ctcomp <- rep(0, length(unique(alls)))
          names(ctcomp) <- unique(alls)
          ctcomp[names(table(alls[c]))] <- table(alls[c])
          ctcomp
        })
      colnames(ctcomp) <- paste0('origin', js.melt[, 2])
      ctcomp.new[rownames(ctcomp), colnames(ctcomp)] <-  ctcomp
    }
    ctcomplist[[pmid]] <- ctcomp.new
  }
  
  saveRDS(ctcomplist, paste0(rdir, 'rm1256_all_simulations_cell_composition.rds'))
}


###############
## test
###############
for (rm.perc in seq(0.1, 0.8, 0.1)){
  print(rm.perc)
  plotdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256_multinomCellComp/', rm.perc, '/plot/')
  rdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256_multinomCellComp/', rm.perc, '/result/')
  
  ctcomplist <- readRDS(paste0(rdir, 'rm1256_all_simulations_cell_composition.rds'))
  
  ## multinome test
  pval <- NULL
  for (ctcomp.new in ctcomplist){
    tryCatch({
      df <- NULL
      for (s in 1:nrow(ctcomp.new)){
        for (b in 1:ncol(ctcomp.new)){
          df <- rbind(df, data.frame(cell = paste0(rownames(ctcomp.new)[s], ';', colnames(ctcomp.new)[b],';', 1:ctcomp.new[s,b]),
                                     branch = rep(colnames(ctcomp.new)[b], ctcomp.new[s,b]),
                                     sample = rep(rownames(ctcomp.new)[s], ctcomp.new[s,b])
          ))
        }
      }
      
      invisible(capture.output(res.bp <- branchPropTest(data = ctcomp.new, design, method = 'multinom')))
      pval = c(pval, res.bp[1,2]) ## first branch
    },error=function(e) {})
  }
  print(summary(pval))
  print(mean(pval < 0.05))
  saveRDS(pval, paste0(rdir, 'rm1256_all_simulations_branch1base_multinom_pvalues2.rds'))
  
  
  ### logit values T test
  # pval <- NULL
  # for (ctcomp in ctcomplist){
  #   ctcomp.logit.tmp <- (ctcomp+1) / (rowSums(ctcomp)+1)
  #   ctcomp <- t(log(ctcomp.logit.tmp/(1-ctcomp.logit.tmp)))
  #   tryCatch({
  #     invisible(capture.output(res.bp <- branchPropTest(data = ctcomp, design, method = 't.test')))
  #     pval = c(pval, res.bp[1]) ## first branch
  #   },error=function(e) {})
  # }
  # summary(pval)
  # saveRDS(pval, paste0(rdir, 'rm1256_all_simulations_logitTtest_pvalues.rds'))
}



###############################################
## summarize reported proportion using pvalues
###############################################
pval <- NULL
for (rm.perc in seq(0.1, 0.8, 0.1)){
  plotdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256_multinomCellComp/', rm.perc, '/plot/')
  rdir <- paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256_multinomCellComp/', rm.perc, '/result/')
  multinom.pval <- readRDS(paste0(rdir, 'rm1256_all_simulations_branch1base_multinom_pvalues2.rds'))
  logitT.pval <- readRDS(paste0(rdir, 'rm1256_all_simulations_logitTtest_pvalues.rds'))
  res <- readRDS(paste0('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/rmBM1256/', rm.perc, '/result/result.rds'))
  names(res)
  str(res)
  pval <- rbind(pval, data.frame(rm.perc = rm.perc * 100, 
                                 multinom.pval = mean(multinom.pval < 0.05),
                                 logitT.pval = mean(logitT.pval < 0.05),
                                 t.pval = res[[4]][1], stringsAsFactors = FALSE))
  
  
  
}
pval <- rbind(pval, c(0,0,0,0))
pval
saveRDS(pval, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/rmBM1256_simulation_branchProp_pvalue.rds')


###########
## plot
##########
pd <- rbind(data.frame(rm.perc = pval[,1], pval = pval[,2], type = 'multinom'),
            data.frame(rm.perc = pval[,1], pval = pval[,3], type = 't-test (logit values)'),
            data.frame(rm.perc = pval[,1], pval = pval[,4], type = 't-test'))
saveRDS(pd, '/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/rmBM1256_simulation_branchProp_pvalue_pd.rds')

str(pd)
pdf('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/rmBM1256_simulation_branchProp_pvalue2.pdf', width = 4, height = 1.9)
ggplot(data = pd, aes(x = rm.perc, y = pval, color = type)) +
  geom_point() + 
  geom_smooth(se=F)+
  theme_classic() + 
  xlab('Reduced cell percentage') + 
  ylab('Report proportion')+
  scale_color_brewer(palette = 'Set2')
dev.off()


pdf('/Users/wenpinhou/Dropbox/trajectory_variability/tree_variability/perf/rmBM1256_simulation_branchProp_pvalue_ttest.pdf', width = 4, height = 1.9)
pd2 <- pd[pd[,'type']!='multinom',]
ggplot(data = pd2, aes(x = rm.perc, y = pval, color = type)) +
  geom_point() + 
  geom_smooth(se=F)+
  theme_classic() + 
  xlab('Reduced cell percentage') + 
  ylab('Report proportion')+
  scale_color_brewer(palette = 'Set1')
dev.off()
