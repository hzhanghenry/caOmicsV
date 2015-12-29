\name{bioNetCircosPlot}
\alias{bioNetCircosPlot}
\title{
    caOmicsV bioNetCircos Layout Plot
}
\description{
    Plot one track of caOmics data on each node of a biological network. 
    Supported plot types include polygon, bar, points, heatmap, and lines.
}
\usage{
    bioNetCircosPlot(dataValues=NULL, plotType="polygon", outer, inner, 
            plotColors=NULL, maxValue=NULL, minValue=NULL)
}
\arguments{
    \item{dataValues}{
    an numeric matrix with columns for samples and rows for genes
    }
    \item{plotType}{
    character vector for plot type, supporting polygon, points, and lines.
    }
    \item{outer}{
    non-negative numeric, outer location of a data track
    }
    \item{inner}{
    non-negative numeric, inner location of a data track
    }
    \item{plotColors}{
    character vector specifying colors for plot items
    }
    \item{maxValue}{
    numeric, the biggest value of plot data or user defined top threshold. Set 
    to NULL to use the biggest value in dataset.
    }
    \item{minValue}{
    numeric, the smallest value of plot data or user defined bottom threshold. 
    Set to NULL to use the smallest value in dataset.
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    data(bionetPlotDemoData)
    expr <- bionetPlotDemoData$heatmapData[[1]]
    bioNet <- bc3net(expr) 
    initializeBioNetCircos(bioNet, totalSamples=60)
    showBioNetNodesLayout()

    binaryData <- bionetPlotDemoData$binaryData[[1]]
    sampleColors <- c(rep("blue", ncol(binaryData)))
    plotType <- "points"
    inner <- 2
    outer <- 3
    bioNetCircosPlot(binaryData, plotType, outer, inner, sampleColors)
}
\keyword{ methods }
