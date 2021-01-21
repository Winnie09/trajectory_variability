GOEnrich <- function(testobj, fdr.cutoff = 0.05,  k = 5, use.clusters = TRUE, type = 'variable', species = 'human'){
  ## k: numeric. number of clusters. Default = 5. Only useful when use.clusters = TRUE. 
  ## type: "time" or "variable". Only useful when "cluster" in not in testobj.
  ## species: currently only work for "human". will include "mouse".
  fdr <- testobj$statistics[,'fdr']
  
  if (sum(fdr < fdr.cutoff) == 0){
    print('There is no differential genes! GoEnrich stopped.')
    break
  } else {
    if (use.clusters){
      if ('cluster' %in% names(testobj)){
        clu <- testobj$cluster    
      } else {
        clu <- clusterGene(testobj, gene = names(testobj$fdr[testobj$fdr < fdr.cutoff]), type = type, k=k)
      }
      diffgeneList <- sapply(unique(clu), function(i){
        names(clu)[clu == i]
      })
    } else (
      diffgeneList <- list(diffgene = rownames(testobj$statistics)[fdr < fdr.cutoff])
    )
    
    resList <- lapply(diffgeneList, function(diffgene){
      allgene <- rownames(testobj$expr.ori)
      if (!'topGO' %in% (.packages())) suppressMessages(library(topGO))
      gl <- diffgene
      back <- allgene
      geneList <- factor(as.integer(back %in% gl))
      names(geneList) <- back
      suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Symbol")
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)})
      sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
      sigres <- sigres[sigres$Annotated >= 10,]
      sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
      ptcount <- 0
      fc <- ((sigres[,"Significant"]+ptcount)/(sum(GOdata@allScores[GOdata@feasible]==1)+ptcount))/((sigres[,"Annotated"]+ptcount)/(sum(GOdata@feasible)+ptcount))
      sigres <- data.frame(sigres,FC=fc)
      sigres <- sigres[order(sigres$FDR,-sigres$FC),]
    })
    return(resList)
  }
    
}


