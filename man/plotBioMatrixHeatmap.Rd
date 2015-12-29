\name{plotBioMatrixHeatmap}
\alias{plotBioMatrixHeatmap}
\title{
    Heatmap Plot on caOmicsV bioMatrix Layout
}
\description{
    Headmap plot method for bioMatrix layout. This function will plot all 
    rows of the input data. bioMatrix layout and graphic device must be 
    initialized first.
}
\usage{
    plotBioMatrixHeatmap(exprData, topAdjust=0, bottomAdjust=0, 
            maxValue=NULL, minValue=NULL, 
            heatmapColor="BlueWhiteRed", skipPlotRow=0)
}
\arguments{
    \item{exprData}{
    numeric matrix (log2 values) with row names
    }
    \item{topAdjust}{
    non-negative numeric, height of top y coordinate should be reduced to show 
    different layers, default 0
    }
    \item{bottomAdjust}{
    non-negative numeric, height of bottom y coordinate should be reduced for 
    a small rectangle, default 0
    }
    \item{maxValue}{
    numeric, value for highest color in heatmap, set to NULL to use the
    maximum value in expression dataset
    }
    \item{minValue}{
    numeric, value for lowest color in heatmap, set to NULL to use the minimum 
    value in expression dataset
    }
    \item{heatmapColor}{
    character vector,one of "BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed", 
    "GreenBlackRed" , or "YellowToRed"
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

    data(biomatrixPlotDemoData)
    exprData <- biomatrixPlotDemoData$heatmapData[[1]]
    exprData <- matrix(exprData[1, 1:50], nrow=1)
    plotBioMatrixHeatmap(exprData, maxValue=3, minValue=-3)
}
\keyword{ methods }
