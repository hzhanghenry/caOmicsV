\name{initializeBioNetCircos}
\alias{initializeBioNetCircos}
\title{
    Initialize caOmicsV bioNetCircos Layout Plot
}
\description{
    Set up parameters for layout of caOmicsV bioNetCircos plot including total 
    number of samples, default width of a sample on a circos track , node 
    radius, padding width between two nodes, plot area of each node, default x 
    and y coordinates of points on a circle with radius of 1, and node layout.
}
\usage{
    initializeBioNetCircos(bioNet, totalSamples=100, sampleWidth=100, 
    nodeRadius=1, nodePadding=1, plotAreaWidth=1, 
    layout=layout.fruchterman.reingold(bioNet))
}
\arguments{
    \item{bioNet}{an igraph object}
    \item{totalSamples}{
    non-negative numeric, total number of samples to be plotted
    }
    \item{sampleWidth}{
    non-negative numeric, total number of points to represent a sample on a 
    circular track.
    }
    \item{nodeRadius}{
    non-negative numeric,  radius of a node on biological network.
    }
    \item{nodePadding}{
    non-negative numeric, padding width between two node on biological 
    network.
    }
    \item{plotAreaWidth}{
    non-negative numeric, outside boundary of plot area of a node, relative 
    to node radius, default 1 (same as nodeRadius)
    }
    \item{layout}{
    two dimentional numeric matrix, node layout of the igraph returned by 
    igraph layout method
    }
}
\author{
    Henry Zhang
}
\value{
    None
}
\examples{
    data(bionetPlotDemoData)
    expr <- bionetPlotDemoData$heatmapData[[1]]
    bioNet <- bc3net(expr)

    initializeBioNetCircos(bioNet, totalSamples=100, sampleWidth=100, 
            nodeRadius=1, nodePadding=1, plotAreaWidth=1, 
            layout=layout.fruchterman.reingold(bioNet)) 
}
\keyword{ methods }
