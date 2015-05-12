\name{getBioMatrixDataRowTop}
\alias{getBioMatrixDataRowTop}
\title{
    Get y Coordinate for Top of A Row on bioMatrix Layout
}
\description{
    Calculate the y coordinate of a row top on bioMatrix layout. The bioMatrix 
    layout must be initialized first
}
\usage{
    getBioMatrixDataRowTop(rowNumber, areaName=c("omicsData", "phenotype"))
}
\arguments{
    \item{rowNumber}{
    non-negative integer, number of the row
    }
    \item{areaName}{
    character vector, either "phenotype" or "omicsdata"
    }
}
\value{
    non-negative numeric, the y coordinate of the row top.
}
\author{
    Henry Zhang
}
\examples{
    initializeBioMatrixPlot()
    yTop <- getBioMatrixDataRowTop(2, areaName="omicsData")
}
\keyword{ methods }
