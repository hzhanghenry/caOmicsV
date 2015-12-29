\name{plotBioNetLines}
\alias{plotBioNetLines}
\title{
    Line Plot on caOmicsV bioNetCircos Layout
}
\description{
    Line plot on caOmicsv bioNetCircos layout. This method plots one track of 
    line for every node. bionetCircos layout and graphic device must be 
    initialized first.
}
\usage{
    plotBioNetLines(dataValues, outer, inner, maxValue=NULL, minValue=NULL, 
            plotColors=rep("black", ncol(dataValues)))
}
\arguments{
    \item{dataValues}{ numeric matrix of plot data }
    \item{outer}{
    non-negative numeric, the outside boundary of plot area from node center
    }
    \item{inner}{
    non-negative numeric, the inside boundary of plot area from node center
    }
    \item{maxValue}{ numeric, the biggest value of plot data }
    \item{minValue}{ numeric, the smallest value of plot data }
    \item{plotColors}{
    character vector or vector of R specification, colours for each sample
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
    plotBioNetLines(expr, outer=13, inner=2)
}
\keyword{ methods }
