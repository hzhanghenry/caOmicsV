\name{plotBioNetBars}
\alias{plotBioNetBars}
\title{
    Bar Plot on caOmicsV bioNetCircos Layout
}
\description{
    Bar plot method for caOmicsV bioNetCircos layout. This will plot one track 
    of bars for every node. bioNetCircos layout and graphic device must be 
    initialized first.
}
\usage{
    plotBioNetBars(dataValues, outer, inner, plotColors)
}
\arguments{
    \item{dataValues}{
    numeric matrix with range of 0 ~ 1 for bar height. total rows of the 
    matrix must be same as the number of nodes and row names must be same as 
    the vertex names in bioNetGraph
    }
    \item{outer}{
    non-negative numeric, the outside boundary of plot area from node center
    }
    \item{inner}{
    non-negative numeric, the inside boundary of plot area from node center
    }
    \item{plotColors}{
    character vector or vector of R color specification, color names for each 
    sample, pre-generated to control sample colors.
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

    methyl <- bionetPlotDemoData$categoryData[[1]]
    sampleColors <- c(rep("red", 20), rep("blue", 20), rep("cyan", 20))
    plotBioNetBars(methyl, outer=1.6, inner=1.5, plotColors=sampleColors)
}
\keyword{ methods }
