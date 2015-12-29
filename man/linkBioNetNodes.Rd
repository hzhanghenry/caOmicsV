\name{linkBioNetNodes}
\alias{linkBioNetNodes}
\title{
    Draw A Customized Arrow between Two Nodes
}
\description{
    Plot a customized arrow between two nodes to replace the edge. 
    bioNetCorcos layout must be initialized first.
}
\usage{
    linkBioNetNodes(fromNode, toNode, lineColor = "black", arrowSize = 1)
}
\arguments{
    \item{fromNode}{ non negative integer, the start node to be linked }
    \item{toNode}{ non negative integer, the end node to be linked }
    \item{lineColor}{ character vector, color of the arrow }
    \item{arrowSize}{
    non-negative numeric, scaling factor for arrow size, default 1
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
    showBioNetNodesLayout()

    linkBioNetNodes(fromNode=2, toNode=5, lineColor="red", arrowSize=1)
}
\keyword{ methods }
