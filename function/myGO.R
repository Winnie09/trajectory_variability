myGO <- function(diffgene, allgene, species = 'human'){
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
  return(sigres)
}
