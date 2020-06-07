plotGene <- function(testptObj, Gene, Mat, Pseudotime, Cellanno, Design,  Alpha=1, Size=0.5, PlotPoints = FALSE, FreeScale = FALSE, BySample = FALSE, type = 'Variable'){
  ## Mat: gene by cell Matrix
  ## testptObj: the output of function testpt() which is a list containing fdr, etc..
 if (type == 'Variable'){
    plotGene_Variable(testptObj=testptObj, Gene=Gene, Mat=Mat, Pseudotime=Pseudotime, Cellanno=Cellanno, Design=Design,  Alpha=Alpha, Size=Size, PlotPoints = PlotPoints, FreeScale = FreeScale, BySample = BySample)
 } else if (type == 'Time'){
    plotGene_Time(testptObj=testptObj, Gene=Gene, Mat=Mat, Pseudotime=Pseudotime, Cellanno=Cellanno, Design=Design,  Alpha=Alpha, Size=Size, PlotPoints = PlotPoints, FreeScale = FreeScale, BySample = BySample)
 }
}
