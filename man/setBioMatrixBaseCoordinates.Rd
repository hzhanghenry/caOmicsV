\name{setBioMatrixBaseCoordinates}
\alias{setBioMatrixBaseCoordinates}
\title{
    Set up Base Coordinates for caOmicsV bioMatrix Layout
}
\description{
    Initialize x and y coordinates for each plot area on bioMatrix layout. 
    This function is for internal use only.
}
\usage{
    setBioMatrixBaseCoordinates(numOfSamples, sampleWidth, columnPadding, 
            sampleHeight, geneNameWidth)
}
\arguments{
    \item{numOfSamples}{
    non-negative integer, number of samples to be plotted
    }
    \item{sampleWidth}{
    non-negative numeric, width of rectangle in inch for each sample, 
    default 0.1
    }
    \item{columnPadding}{
    non-negative numeric, padding width in inch between two samples, 
    default 0.025
    }
    \item{sampleHeight}{
    non-negative numeric, height of rectangle in inchfor each sample, 
    default 0.4
    }
    \item{geneNameWidth}{
    non-negative numeric, width of gene labeling area in inch, default 1
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    #  Internal use only
}
\keyword{ methods }
