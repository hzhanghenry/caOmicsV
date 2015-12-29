\name{showBioNetNodesLayout}
\alias{showBioNetNodesLayout}
\title{
    Display Nodes Layout of caOmicsV bioNetCircos Plot
}
\description{
    Display bioNetCircos layout with empty nodes and edges in order to 
    optimize bioNetCircos layout. An igraph object must be initialized first.
}
\usage{
    showBioNetNodesLayout(bgColor=grey(0.75, alpha=0.5))
}
\arguments{
    \item{bgColor}{
    character vector for a color name or a R color specification 
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
    showBioNetNodesLayout(bgColor=grey(0.75, alpha=0.5))
}
\keyword{ methods }
