\name{linkBioNetSamples}
\alias{linkBioNetSamples}
\title{
    Link Two samples Inside a Node with Quadratic Bezier Curve
}
\description{
    Draw a quadratic Bezier curve line between two samples inside of a node. 
    bioNetCorcos layout must be initialized first.
}
\usage{
    linkBioNetSamples(nodeIndex, fromSample, toSample, outer, plotColors)
}
\arguments{
    \item{nodeIndex}{
    non-negative integer, the node on which link line is drawn
    }
    \item{fromSample}{
    non-negative integer, the first sample to be linked
    }
    \item{toSample}{
    non-negative integer, the second sample to be linked
    }
    \item{outer}{
    non-negative numeric, the start and end of link line relative to node 
    center
    }
    \item{plotColors}{
    character vector, color for the link line (ribbon)
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

    linkBioNetSamples(nodeIndex=2, fromSample=10, toSample=20, outer=0.9, 
            plotColors="red")
}
\keyword{ methods }
