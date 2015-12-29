\name{getBioNetPlotLocations}
\alias{getBioNetPlotLocations}
\title{
    Get bioNetCircos Plot Locations
}
\description{
    Get plot locations for a node on caOmicsC bioNetCircos layout
}
\usage{
    getBioNetPlotLocations(nodeCenter, outer, inner)
}
\arguments{
    \item{nodeCenter}{
    numeric, x and y coordinates of the node center
    }
    \item{outer}{
    non-negative numeric, outer limit of plot track relative to node center
    }
    \item{inner}{
    non-negative numeric, inner limit of plot track relative to ode center
    }
}
\value{
    \item{ nodeCenter }{ numeric, x and y coordinates of the node center }
    \item{ outPositions }{ 
    two dimensional numeric matrix for x and y coordinates of outer boundary 
    for plot 
    }
    \item{ inPositions }{ 
    two dimensional numeric matrix for x and y coordinates of inner boundary 
    for plot 
    }
    \item{ positionIndex }{ 
    matrix with index of x and y coordinates of base plot position for 
    each sample
    }
}
\author{
    Henry Zhang
}
\examples{
    data(bionetPlotDemoData)
    expr <- bionetPlotDemoData$heatmapData[[1]]
    bioNet <- bc3net(expr)

    initializeBioNetCircos(bioNet)

    bioNetGraph <- getBioNetGraph()
    outer <- 1.5
    inner <- 1.0

    nodeCenter <- as.numeric(bioNetGraph$layout[1,])
    plotLocations <- getBioNetPlotLocations(nodeCenter, outer, inner)
}
\keyword{ methods }
