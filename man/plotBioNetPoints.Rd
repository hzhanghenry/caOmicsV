\name{plotBioNetPoints}
\alias{plotBioNetPoints}
\title{
    Point Plot on caOmicsV bioNetCircos Layout
}
\description{
    Point plot on caOmicsv bioNetCircos layout. This method plots one track of 
    points for every node. bioNetCircos layout and graphic device must be 
    initialized first.
}
\usage{
    plotBioNetPoints(dataValues, maxValue=NULL, 
            minValue=NULL, outer, inner, 
            plotColors=rep("black", ncol(dataValues)), 
            sizeByValue=FALSE, pch=".")
}
\arguments{
    \item{dataValues}{ numeric matrix of plot data }
    \item{maxValue}{ numeric, the biggest value of plot data }
    \item{minValue}{ numeric, the smallest value of plot data }
    \item{outer}{ 
    non-negative numeric, the outer boundary of plot area from node center
    }
    \item{inner}{
    non-negative numeric, the inner boundary of plot area from node center
    }
    \item{plotColors}{
    character vector or vector of R color specification, colors for each 
    sample
    }
    \item{sizeByValue}{
    logic, if true, the data value will  be used for point size
    }
    \item{pch}{ character for point type, same as the one in par() }
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
    plotBioNetPoints(expr, outer=3, inner=2)
}
\keyword{ methods }
