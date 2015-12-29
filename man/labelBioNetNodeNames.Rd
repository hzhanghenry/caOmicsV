\name{labelBioNetNodeNames}
\alias{labelBioNetNodeNames}
\title{
    Label Names for Each Node on Network Graph
}
\description{
    Plot name for one or more nodes on an igraph network. igraph object and 
    graphic device must be initialized first.
}
\usage{
    labelBioNetNodeNames(nodeList=NULL, labelColor="black", 
        labelLocation=c("bottom", "left", "top", "right"), 
        labelOffset=0.5)
}
\arguments{
    \item{nodeList}{
    non-negative integer, index of node(s) on a bioNet layout
    }
    \item{labelColor}{ character vector, colored for text (labels) }
    \item{labelLocation}{
    character vector, location relative to node center, either "bottom", 
    "left", "top", or "right"
    }
    \item{labelOffset}{
    non-negative numeric, distance from node outside boundary
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
    labelBioNetNodeNames(nodeList=c(1, 4, 7), labelColor="black", 
        labelLocation="bottom")
}
\keyword{ methods }
