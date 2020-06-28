plotGene <- function(testptObj, Gene, Mat, Pseudotime, Cellanno, Design,  Alpha=1, Size=0.5, PlotPoints = FALSE, FreeScale = FALSE, BySample = FALSE, type = 'Variable'){
  ## Mat: gene by cell Matrix
  ## testptObj: the output of function testpt() which is a list containing fdr, etc..
  ### pseudotime: a dataframe of 1st column cell, 2nd column pseudotime
  if (ncol(Pseudotime) > 1) {
    Pseudotime <- data.frame(Cell = Pseudotime[,1], Pseudotime = as.numeric(Pseudotime[,2]), stringsAsFactors = FALSE)
    Pseudotime <- Pseudotime[order(Pseudotime[,2]), ]
    psn <- Pseudotime[,2]
    names(psn) <- Pseudotime[,1]
    Pseudotime <- psn
  } else {
    print('Note: pseudotime should be a dataframe containing 1st column cell and 2nd column pseudotime.')
  } 
 if (type == 'Variable'){
    plotGene_Variable(testptObj=testptObj, Gene=Gene, Mat=Mat, Pseudotime=Pseudotime, Cellanno=Cellanno, Design=Design,  Alpha=Alpha, Size=Size, PlotPoints = PlotPoints, FreeScale = FreeScale, BySample = BySample)
 } else if (type == 'Time'){
    plotGene_Time(testptObj=testptObj, Gene=Gene, Mat=Mat, Pseudotime=Pseudotime, Cellanno=Cellanno, Design=Design,  Alpha=Alpha, Size=Size, PlotPoints = PlotPoints, FreeScale = FreeScale, BySample = BySample)
 }
}
