\name{getBioNetSamplePlotPosition}
\alias{getBioNetSamplePlotPosition}
\title{
    Calculate x and y Coordinates for Each Sample on Default Node
}
\description{
    Calculate x and y coordinates for each sample on default node. The output 
    will be a three column matrix representing the left, center, and right 
    position for each sample on circumference of default node. The center 
    positions are for points plot and others are for polygon plot. The
    bionetCircos layout must be initialized first.
}
\usage{
    getBioNetSamplePlotPosition(totalSamples)
}
\arguments{
    \item{totalSamples}{
    non-negative integer, total numbe of samples to be plotted
    }
}
\value{
    Matrix with index of x and y coordinates for each sample
}
\author{
    Henry Zhang
}
\examples{
    totalSamples <- 100
    samplePositions <- getBioNetSamplePlotPosition(totalSamples)
}
\keyword{ methods }
