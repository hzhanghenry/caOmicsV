\name{getBioNetNodeLinkLine}
\alias{getBioNetNodeLinkLine}
\title{
    Get X and Y Coordinates for An Arrow between Two Nodes
}
\description{
    Calculate x and y coordinates for an customized arrow head and tail with 
    defined length to connect two nodes.
}
\usage{
    getBioNetNodeLinkLine(lineX, lineY, arrowSize=1, lineLength)
}
\arguments{
    \item{lineX}{
    numeric vector, x coordinates of the link line
    }
    \item{lineY}{
    numeric vector, y coordinates of the link line
    }
    \item{arrowSize}{
    non-negative numeric, scaling factor for arrow size, default 1
    }
    \item{lineLength}{
    non negative integer, the length of link line
    }
}
\details{
    An arrow is drawn as a polygon. By default, the arrow is in inside of a 
    circle (radius 1) without tail and it points to radian 0. The tail, if 
    any, will be added to the left.
}
\value{
    A two dimensional numeric matrix for x and y coordinates of the arrow.
}
\author{
    Henry Zhang
}
\examples{
    from <- c(1, 1)
    to   <- c(2, 2)
    lineX <- seq(from[1], to[1], 1000)
    lineY <- seq(from[2], to[2], 1000)

    lineLength <- sqrt((from[1]-to[1])^2 + (from[2]-to[2])^2)
    positions <- getBioNetNodeLinkLine(lineX, lineY, arrowSize=1, lineLength)
}
\keyword{ methods }
