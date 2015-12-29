\name{plotBioNetHeatmap}
\alias{plotBioNetHeatmap}
\title{
    Heatmap Plot for Each Node on caOmicsV bioNetCircos Layout
}
\description{
    Headmap plot on caOmicsv bioNetCircos layout. This method plots one track 
    of heatmap for every node. bioNetCircos layout and graphic device must be 
    initialized first.
}
\usage{
    plotBioNetHeatmap(dataValues, maxValue=NULL, minValue=NULL, 
            outer, inner, plotColors)
}
\arguments{
    \item{dataValues}{
    numeric matrix of log2 values for heatmap plot. Total rows of the matrix 
    must be same as the number of nodes and rownames must be same as the 
    vertex names in bioNetGraph
    }
    \item{maxValue}{ numeric, the biggest value of plot data }
    \item{minValue}{ numeric, the smallest  value of plot data }
    \item{outer}{
    non-negative numeric, the outside boundary of plot area from node center
    }
    \item{inner}{
    non-negative numeric, the inside boundary of plot area close to node 
    center
    }
    \item{plotColors}{
    character vector, one of "BlueWhiteRed", "GreenWhiteRed", 
    "GreenYellowRed", "GreenBlackRed", "YellowToRed", "BlackOnly".
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
    plotBioNetHeatmap(expr, outer=3, inner=2, plotColors="BlueWhiteRed")
}
\keyword{ methods }
