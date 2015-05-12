\name{plotBioMatrixSampleNames}
\alias{plotBioMatrixSampleNames}
\title{
    Label Sample Names on the Top of caOmicsV bioMatrix Layout
}
\description{
    Text plot on the top of phenotype area to show sample names. bioMatrix 
    layout and graphic device must be initialized first.
}
\usage{
    plotBioMatrixSampleNames(sampleNames, sampleColors)
}
\arguments{
    \item{sampleNames}{
    character vector, sample names to be plotted
    }
    \item{sampleColors}{
    character vector or R color name(s) for text color(s)
    }
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

    sampleNames <- paste("Sample", 1:50)
    sampleColors <- c(rep("red", 30), rep("blue", 30))

    plotBioMatrixSampleData(1, areaName="phenotype", fillColor = sampleColors)
    plotBioMatrixSampleNames(sampleNames, sampleColors)
}
\keyword{ methods }
