\name{setBioNetNodeLayout}
\alias{setBioNetNodeLayout}
\title{
    Set up Layout for the igraph Object
}
\description{
    Set up node layout on a biological network. The layout is taken from an 
    igraph layout and scaled to allocate circos plot area for each node.
}
\usage{
    setBioNetNodeLayout(bioNet, layout=layout.auto(bioNet))
}
\arguments{
    \item{bioNet}{
    an igraph object for which layout will be set
    }
    \item{layout}{
    layout object from layout method of igraph package
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
    setBioNetNodeLayout(bioNet, layout=layout.auto(bioNet)) 
}
\keyword{ methods }
