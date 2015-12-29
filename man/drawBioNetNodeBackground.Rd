\name{drawBioNetNodeBackground}
\alias{drawBioNetNodeBackground}
\title{
    Draw Background for A Data Track On Nodes of caOmicsV bioNet Layout
}
\description{
    Paint (with any color other than white) background for a circular track 
    on each node on caOmicsV bioNetCircos layout. Graphic device and igraph 
    object must exist.
}
\usage{
    drawBioNetNodeBackground(trackLocations, bgColor=gray(0.9, alpha=0.5))
}
\arguments{
    \item{trackLocations}{
    a list returned by getBioNetPlotLocations() function.
    }
    \item{bgColor}{
    character vector or R color specification for background color
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
    initializeBioNetCircos(bioNet)

    bioNetGraph <- getBioNetGraph()
    outer <- 1.5
    inner <- 1.0

    nodeCenter <- as.numeric(bioNetGraph$layout[1,])
    plotLocations <- getBioNetPlotLocations(nodeCenter, outer, inner)

    showBioNetNodesLayout()
    eraseBioNetNode()
    drawBioNetNodeBackground(plotLocations)
}
\keyword{ methods }

