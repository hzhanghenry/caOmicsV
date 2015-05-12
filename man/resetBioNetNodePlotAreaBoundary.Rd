\name{resetBioNetNodePlotAreaBoundary}
\alias{resetBioNetNodePlotAreaBoundary}
\title{
    Update Node Plot Area Boundary on caOmicsV bioNetCircos Layout
}
\description{
    Record the plotted area boundary for all nodes on the igraph. These 
    boundary may be needed for drawn customized arrows and label node names. 
    bioNetCircos layout must be initialized first. This function is for 
    internal use only.
}
\usage{
    resetBioNetNodePlotAreaBoundary(inner=getBioNetNodeRadius(), outer)
}
\arguments{
    \item{inner}{
    non-negative numeric, the inner boundary of area that has been plotted
    }
    \item{outer}{
    non-negative numeric, the outer boundary of area that has been plotted
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    # Internal use only
}
\keyword{ methods }
