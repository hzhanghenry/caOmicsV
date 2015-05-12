\name{setBioMatrixPlotArea}
\alias{setBioMatrixPlotArea}
\title{
    Set up Plot Area for caOimcsV bioMatrix Layout
}
\description{
    Set up plot area including of sample name area (sampleHeight), phenotype 
    area, gene label area (geneNameWidth), and remark area (remarkWidth, and 
    other descriptions). bioMatrix layout muse be initialized first. 
}
\usage{
    setBioMatrixPlotArea()
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    initializeBioMatrixPlot(numOfGenes=1, numOfSamples=50)
    setBioMatrixPlotArea()
}
\keyword{ methods }
