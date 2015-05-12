\name{setBioNetPlotParameters}
\alias{setBioNetPlotParameters}
\title{
    Initialize Plot Parameters for caOmicsV bioNetCircos Plot
}
\description{
    Set up plot parameters for bioNetCircos layout including of totalSamples, 
    sampleWidth, nodeRadius, nodePadding, plotAreaWidth, outer and inner 
    boundaries of plotArea, as well as default plot colors. This function is 
    for internal used only.
}
\usage{
    setBioNetPlotParameters(totalSamples, sampleWidth, nodeRadius, 
                    nodePadding, plotAreaWidth)
}
\arguments{
    \item{totalSamples}{
    non-negative integer, total number of samples to be plotted
    }
    \item{sampleWidth}{
    non-negative integer, number of points of a sample on the circular line
    }
    \item{nodeRadius}{
    non-negative numeric, radius of the node
    }
    \item{nodePadding}{
    non-negative numeric, padding width in inch between two nodes
    }
    \item{plotAreaWidth}{
    non-negative numeric, total width of plot area.
    }
}
\value{
    None.
}
\author{
    Henry zhang
}
\examples{
    #  Internal used only.
}
\keyword{ methods }
