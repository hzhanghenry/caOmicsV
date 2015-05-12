\name{plotBioMatrixBars}
\alias{plotBioMatrixBars}
\title{
    Bar Plot on caOmicsV bioMatrix Layout}
\description{
    Bar plot method for caOmicsV bioMatrix layout with non-negative numeric 
    matrix or vector with values in range of 0 ~ 1. bioMatrix layout and 
    graphic device must be initialized first.
}
\usage{
    plotBioMatrixBars(barData, barColor="red", areaName="omicsData", 
    byRow=TRUE, skipPlotRow=0, skipPlotColumns=0)
}
\arguments{
    \item{barData}{
    non-negative numeric matrix or vector with values in range of 0 ~ 1
    }
    \item{barColor}{
    character vector for color name or R color specification
    }
    \item{areaName}{
    character vector, name of plot area, currentlt use "omicsData" only
    }
    \item{byRow}{ logic, whether plot bars for each row or not }
    \item{skipPlotRow}{
    non-negative integer, how many row(s) to be skipped from first row
    }
    \item{skipPlotColumns}{
    non-negative integer, how many row(s) to be skipped from first column
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

    barData <- matrix(c(rep(0.25, 15), rep(0.75, 20), rep(0.5, 15)), nrow=1)
    plotBioMatrixBars(barData, barColor="red") 
}
\keyword{ methods }
