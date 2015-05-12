\name{eraseBioNetNode}
\alias{eraseBioNetNode}
\title{
    Erase Background of All Nodes on caOmicsV bioNetCircos Layout
}
\description{
    Erase all nodes on a caOmicsV bioNetCircos layout except of the edges. 
    Graphic device and igraph object must be initialized first.
}
\usage{
    eraseBioNetNode()
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    data(caOmicsV.bionet.eset)
    expr <- caOmicsV.bionet.eset$heatmapData[[1]]
    bioNet <- bc3net(expr)
    initializeBioNetCircos(bioNet)

    showBioNetNodesLayout()
    eraseBioNetNode()
}
\keyword{ methods }
