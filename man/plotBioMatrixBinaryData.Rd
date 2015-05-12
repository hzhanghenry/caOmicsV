\name{plotBioMatrixBinaryData}
\alias{plotBioMatrixBinaryData}
\title{
    Binary Data Plot on caOmicsV bioMatrix Layout
}
\description{
    Plot binary data as points in the inside of each rectangle(sample). This 
    function plot all rows on omics data area and only the positive samples 
    will be shown with colored points. For one row plot, pass data as vector 
    and supply correct skipPlotRow parameter to define where to plot. bioMatrix 
    layout and graphic device must be initialized first.
}
\usage{
    plotBioMatrixBinaryData(binaryData, areaName="omicsData", scatterType=19, 
        scatterSize=1, totalSubRow=1, subRowIndex=1, sampleColor="black", 
        skipPlotRow=0)
}
\arguments{
    \item{binaryData}{
    vector or matrix with values of 0 and 1 only 
    }
    \item{areaName}{ character vector, either "omicsData" or "phenotype" }
    \item{scatterType}{ non-negative integer, same as pch, default 19 }
    \item{scatterSize}{ non-negative numeric, same as cex }
    \item{totalSubRow}{ 
    non-negative integer, how many sub-rows in a sample area 
    }
    \item{subRowIndex}{ non-negative integer, which subrow will be plotted }
    \item{sampleColor}{
    character vector for color name(s) or R color specification 
    }
    \item{skipPlotRow}{
    non-negative integer, total rows on plot area that should be skipped, 
    default 0
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
    showBioMatrixPlotLayout("Gene", paste("Sample", 1:50),  "Diagnosis")

    binaryData <- matrix(c(rep(1, 15), rep(0, 20), rep(1, 15)), nrow=1)
    plotBioMatrixBinaryData(binaryData, scatterType=16) 
}
\keyword{ methods }
